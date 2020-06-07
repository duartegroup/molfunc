from functools import wraps
from molfunc.exceptions import NoAtomsInMolecule


def requires_atoms():
    """A function requiring atoms to run"""

    def func_decorator(func):
        @wraps(func)
        def wrapped_function(*args, **kwargs):
            molecule = args[0]

            assert hasattr(args[0], 'n_atoms')
            assert hasattr(args[0], 'atoms')

            if molecule.atoms is None or molecule.n_atoms == 0:
                raise NoAtomsInMolecule

            return func(*args, **kwargs)

        return wrapped_function
    return func_decorator
