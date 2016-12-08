class Class_Operation(object):
    """
    Class Operation
    """

    def __init__(self, operation, cls1, cls2):
        """
        """
        self.cls1 = cls1
        self.cls2 = cls2
        self._operation = operation

    def __getattr__(self, attr):
        try:
            return getattr(getattr(self.cls1, attr),
                           self._operation)(getattr(self.cls2, attr))
        except AttributeError:
            return lambda *args: getattr(
                getattr(self.cls1, attr)(*args),
                self._operation)(getattr(self.cls2, attr)(*args))

    def __add__(self, other):
        return Class_Operation('__add__', self, other)

    def __sub__(self, other):
        return Class_Operation('__sub__', self, other)

    def __mul__(self, other):
        return Class_Operation('__mul__', self, other)

    def __gt__(self, other):
        return Class_Operation('__gt__', self, other)
