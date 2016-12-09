"Some utility methods"


class ClassOperationResult(object):
    """
    Result of Class Operation
    """
    def __init__(self, operation, cls1, cls2):
        """
        """
        self.cls1 = cls1
        self.cls2 = cls2
        self._operation = operation

    def __getattr__(self, attr):
        attr1 = getattr(self.cls1, attr)
        attr2 = getattr(self.cls2, attr)

        try:
            return getattr(attr1, self._operation)(attr2)

        except AttributeError:
            return (lambda *args, **kwargs:
                    getattr(attr1(*args), self._operation)(attr2(*args)))


class ClassOperation(object):
    """
    Define arithmetic operation on Class
    """
    def __add__(self, other):
        return ClassOperationResult('__add__', self, other)

    def __sub__(self, other):
        return ClassOperationResult('__sub__', self, other)

    def __mul__(self, other):
        return ClassOperationResult('__mul__', self, other)

    def __gt__(self, other):
        return ClassOperationResult('__gt__', self, other)
