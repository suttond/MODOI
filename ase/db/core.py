import collections
import functools
import operator
import os
import re
from time import time

from ase.atoms import Atoms, symbols2numbers
from ase.db.row import atoms2dict, AtomsRow
from ase.calculators.calculator import all_properties, all_changes
from ase.data import atomic_numbers
from ase.parallel import world, broadcast, DummyMPI
from ase.utils import Lock, basestring


T2000 = 946681200.0  # January 1. 2000
YEAR = 31557600.0  # 365.25 days


def now():
    """Return time since January 1. 2000 in years."""
    return (time() - T2000) / YEAR
        

seconds = {'s': 1,
           'm': 60,
           'h': 3600,
           'd': 86400,
           'w': 604800,
           'M': 2629800,
           'y': YEAR}

longwords = {'s': 'second',
             'm': 'minute',
             'h': 'hour',
             'd': 'day',
             'w': 'week',
             'M': 'month',
             'y': 'year'}

ops = {'<': operator.lt,
       '<=': operator.le,
       '=': operator.eq,
       '>=': operator.ge,
       '>': operator.gt,
       '!=': operator.ne}

invop = {'<': '>=', '<=': '>', '>=': '<', '>': '<=', '=': '!=', '!=': '='}

word = re.compile('[_a-zA-Z][_0-9a-zA-Z]*$')

reserved_keys = set(all_properties + all_changes +
                    ['id', 'unique_id', 'ctime', 'mtime', 'user',
                     'momenta', 'constraints',
                     'calculator', 'calculator_parameters',
                     'key_value_pairs', 'data'])

numeric_keys = set(['id', 'energy', 'magmom', 'charge', 'natoms'])


def check(key_value_pairs):
    for key, value in key_value_pairs.items():
        if not word.match(key) or key in reserved_keys:
            raise ValueError('Bad key: {0}'.format(key))
        if not isinstance(value, (int, float, basestring)):
            raise ValueError('Bad value: {0}'.format(value))

            
def connect(name, type='extract_from_name', create_indices=True,
            use_lock_file=True, append=True):
    """Create connection to database.
    
    name: str
        Filename or address of database.
    type: str
        One of 'json', 'db', 'postgresql', 'mysql'
        (JSON, SQLite, PostgreSQL, MySQL/MariaDB).
        Default is 'extract_from_name', which will ... guess the type
        from the name.
    use_lock_file: bool
        You can turn this off if you know what you are doing ...
    append: bool
        Use append=False to start a new database.
    """
    
    if type == 'extract_from_name':
        if name is None:
            type = None
        elif name.startswith('pg://'):
            type = 'postgresql'
        else:
            type = os.path.splitext(name)[1][1:]

    if type is None:
        return Database()

    if not append and world.rank == 0 and os.path.isfile(name):
        os.remove(name)
        
    if type == 'json':
        from ase.db.jsondb import JSONDatabase
        return JSONDatabase(name, use_lock_file=use_lock_file)
    if type == 'db':
        from ase.db.sqlite import SQLite3Database
        return SQLite3Database(name, create_indices, use_lock_file)
    if type == 'postgresql':
        from ase.db.postgresql import PostgreSQLDatabase
        return PostgreSQLDatabase(name[5:])
    raise ValueError('Unknown database type: ' + type)


def lock(method):
    """Decorator for using a lock-file."""
    @functools.wraps(method)
    def new_method(self, *args, **kwargs):
        if self.lock is None:
            return method(self, *args, **kwargs)
        else:
            with self.lock:
                return method(self, *args, **kwargs)
    return new_method


def parallel(method):
    """Decorator for broadcasting from master to slaves using MPI."""
    if world.size == 1:
        return method
        
    @functools.wraps(method)
    def new_method(*args, **kwargs):
        ex = None
        result = None
        if world.rank == 0:
            try:
                result = method(*args, **kwargs)
            except Exception as ex:
                pass
        ex, result = broadcast((ex, result))
        if ex is not None:
            raise ex
        return result
    return new_method


def parallel_generator(generator):
    """Decorator for broadcasting yields from master to slaves using MPI."""
    if world.size == 1:
        return generator
        
    @functools.wraps(generator)
    def new_generator(*args, **kwargs):
        if world.rank == 0:
            for result in generator(*args, **kwargs):
                result = broadcast(result)
                yield result
            broadcast(None)
        else:
            result = broadcast(None)
            while result is not None:
                yield result
                result = broadcast(None)
    return new_generator

    
def convert_str_to_float_or_str(value):
    """Safe eval()"""
    try:
        value = float(value)
    except ValueError:
        value = {'True': 1.0, 'False': 0.0}.get(value, value)
    return value
    

class Database:
    """Base class for all databases."""
    def __init__(self, filename=None, create_indices=True,
                 use_lock_file=False):
        if isinstance(filename, str):
            filename = os.path.expanduser(filename)
        self.filename = filename
        self.create_indices = create_indices
        if use_lock_file and isinstance(filename, str):
            self.lock = Lock(filename + '.lock', world=DummyMPI())
        else:
            self.lock = None
            
    @parallel
    @lock
    def write(self, atoms, key_value_pairs={}, data={}, **kwargs):
        """Write atoms to database with key-value pairs.
        
        atoms: Atoms object
            Write atomic numbers, positions, unit cell and boundary
            conditions.  If a calculator is attached, write also already
            calculated properties such as the energy and forces.
        key_value_pairs: dict
            Dictionary of key-value pairs.  Values must be strings or numbers.
        data: dict
            Extra stuff (not for searching).
            
        Key-value pairs can also be set using keyword arguments::
            
            connection.write(atoms, name='ABC', frequency=42.0)
            
        Returns integer id of the new row.
        """
        
        if atoms is None:
            atoms = Atoms()
        
        kvp = dict(key_value_pairs)  # modify a copy
        kvp.update(kwargs)
        
        id = self._write(atoms, kvp, data)
        return id
        
    def _write(self, atoms, key_value_pairs, data):
        check(key_value_pairs)
        return 1

    @parallel
    @lock
    def reserve(self, **key_value_pairs):
        """Write empty row if not already present.
        
        Usage::
            
            id = conn.reserve(key1=value1, key2=value2, ...)
        
        Write an empty row with the given key-value pairs and
        return the integer id.  If such a row already exists, don't write
        anything and return None.
        """
        
        for dct in self._select([],
                                [(key, '=', value)
                                 for key, value in key_value_pairs.items()]):
            return None

        atoms = Atoms()
        
        calc_name = key_value_pairs.pop('calculator', None)
        
        if calc_name:
            # Allow use of calculator key
            assert calc_name.lower() == calc_name
            
            # Fake calculator class:
            class Fake:
                name = calc_name
                
                def todict(self):
                    return {}
                
                def check_state(self, atoms):
                    return ['positions']
            
            atoms.calc = Fake()
            
        id = self._write(atoms, key_value_pairs, {})
        
        return id
        
    def __delitem__(self, id):
        self.delete([id])
        
    def get_atoms(self, selection=None, attach_calculator=False,
                  add_additional_information=False, **kwargs):
        """Get Atoms object.
        
        selection: int, str or list
            See the select() method.
        attach_calculator: bool
            Attach calculator object to Atoms object (default value is
            False).
        add_additional_information: bool
            Put key-value pairs and data into Atoms.info dictionary.
        
        In addition, one can use keyword arguments to select specific
        key-value pairs.
        """
            
        row = self.get(selection, **kwargs)
        return row.toatoms(attach_calculator, add_additional_information)

    def __getitem__(self, selection):
        return self.get(selection)

    def get(self, selection=None, **kwargs):
        """Select a single row and return it as a dictionary.
        
        selection: int, str or list
            See the select() method.
        fancy: bool
            return fancy dictionary with keys as attributes (this is the
            default).
        """

        rows = list(self.select(selection, limit=2, **kwargs))
        if not rows:
            raise KeyError('no match')
        assert len(rows) == 1, 'more than one row matched'
        return rows[0]

    def parse_selection(self, selection, **kwargs):
        if selection is None or selection == '':
            expressions = []
        elif isinstance(selection, int):
            expressions = [('id', '=', selection)]
        elif isinstance(selection, list):
            expressions = selection
        else:
            expressions = [w.strip() for w in selection.split(',')]
        keys = []
        comparisons = []
        for expression in expressions:
            if isinstance(expression, (list, tuple)):
                comparisons.append(expression)
                continue
            if expression.count('<') == 2:
                value, expression = expression.split('<', 1)
                if expression[0] == '=':
                    op = '>='
                    expression = expression[1:]
                else:
                    op = '>'
                key = expression.split('<', 1)[0]
                comparisons.append((key, op, value))
            for op in ['!=', '<=', '>=', '<', '>', '=']:
                if op in expression:
                    break
            else:
                if expression in atomic_numbers:
                    comparisons.append((expression, '>', 0))
                else:
                    keys.append(expression)
                continue
            key, value = expression.split(op)
            comparisons.append((key, op, value))

        cmps = []
        for key, value in kwargs.items():
            comparisons.append((key, '=', value))
            
        for key, op, value in comparisons:
            if key == 'age':
                key = 'ctime'
                op = invop[op]
                value = now() - time_string_to_float(value)
            elif key == 'formula':
                assert op == '='
                numbers = symbols2numbers(value)
                count = collections.defaultdict(int)
                for Z in numbers:
                    count[Z] += 1
                cmps.extend((Z, '=', count[Z]) for Z in count)
                key = 'natoms'
                value = len(numbers)
            elif key in atomic_numbers:
                key = atomic_numbers[key]
                value = int(value)
            elif isinstance(value, basestring):
                value = convert_str_to_float_or_str(value)
            if key in numeric_keys and not isinstance(value, (int, float)):
                msg = 'Wrong type for "{0}{1}{2}" - must be a number'
                raise ValueError(msg.format(key, op, value))
            cmps.append((key, op, value))
            
        return keys, cmps

    @parallel_generator
    def select(self, selection=None, filter=None, explain=False,
               verbosity=1, limit=None, offset=0, sort=None, **kwargs):
        """Select rows.
        
        Return AtomsRow iterator with results.  Selection is done
        using key-value pairs and the special keys:
            
            formula, age, user, calculator, natoms, energy, magmom
            and/or charge.
        
        selection: int, str or list
            Can be:
            
            * an integer id
            * a string like 'key=value', where '=' can also be one of
              '<=', '<', '>', '>=' or '!='.
            * a string like 'key'
            * comma separated strings like 'key1<value1,key2=value2,key'
            * list of strings or tuples: [('charge', '=', 1)].
        filter: function
            A function that takes as input a row and returns True or False.
        explain: bool
            Explain query plan.
        verbosity: int
            Possible values: 0, 1 or 2.
        limit: int or None
            Limit selection.
        """
        
        if sort:
            if sort == 'age':
                sort = '-ctime'
            elif sort == '-age':
                sort = 'ctime'
            elif sort.lstrip('-') == 'user':
                sort += 'name'
            
        keys, cmps = self.parse_selection(selection, **kwargs)
        for row in self._select(keys, cmps, explain=explain,
                                verbosity=verbosity,
                                limit=limit, offset=offset, sort=sort):
            if filter is None or filter(row):
                yield row
                
    def count(self, selection=None, **kwargs):
        n = 0
        for row in self.select(selection, **kwargs):
            n += 1
        return n
        
    @parallel
    @lock
    def update(self, ids, delete_keys=[], block_size=1000,
               **add_key_value_pairs):
        """Update row(s).
        
        ids: int or list of int
            ID's of rows to update.
        delete_keys: list of str
            Keys to remove.
            
        Use keyword argumnts to add new keys-value pairs.
            
        Returns number of key-value pairs added and removed.
        """
        check(add_key_value_pairs)

        if isinstance(ids, int):
            ids = [ids]
            
        B = block_size
        nblocks = (len(ids) - 1) // B + 1
        M = 0
        N = 0
        for b in range(nblocks):
            m, n = self._update(ids[b * B:(b + 1) * B], delete_keys,
                                add_key_value_pairs)
            M += m
            N += n
        return M, N

    def delete(self, ids):
        """Delete rows."""
        raise NotImplementedError

        
def time_string_to_float(s):
    if isinstance(s, (float, int)):
        return s
    s = s.replace(' ', '')
    if '+' in s:
        return sum(time_string_to_float(x) for x in s.split('+'))
    if s[-2].isalpha() and s[-1] == 's':
        s = s[:-1]
    i = 1
    while s[i].isdigit():
        i += 1
    return seconds[s[i:]] * int(s[:i]) / YEAR


def float_to_time_string(t, long=False):
    t *= YEAR
    for s in 'yMwdhms':
        x = t / seconds[s]
        if x > 5:
            break
    if long:
        return '{0:.3f} {1}s'.format(x, longwords[s])
    else:
        return '{0:.0f}{1}'.format(round(x), s)
