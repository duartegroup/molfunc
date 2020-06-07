class XYZfileMalformatted(Exception):
    pass


class XYZfileDidNotExist(Exception):
    pass


class NotXYZfile(Exception):
    pass


class NoAtomsToPrint(Exception):
    pass


class RAtomNotFound(Exception):
    pass


class RDKitFailed(Exception):
    pass


class NoAtomsInMolecule(Exception):
    pass


class MolFuncCritical(Exception):
    def __init__(self, message=''):
        super().__init__(message)


class DatomsNotValid(MolFuncCritical):
    pass


class RAtomInvalidValence(Exception):
    pass


class CombinationFailed(MolFuncCritical):
    pass
