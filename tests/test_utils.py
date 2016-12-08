from ne2001.utils import Class_Operation


class Foo(Class_Operation):
    """
    """

    def __init__(self, func, x):
        self.x = x
        self.func = func

def test():
    f1 = Foo(lambda x: 2*x, 1)
    f2 = Foo(lambda x: x**2, 2)

    f = f1 + f2
    assert f.x == f1.x + f2.x
    assert f.func(3) == f1.func(3) + f2.func(3)

    f = f1 - f2
    assert f.x == f1.x - f2.x
    assert f.func(3) == f1.func(3) - f2.func(3)

    f = f1 + f2
    assert f.x == f1.x + f2.x
    assert f.func(3) == f1.func(3) + f2.func(3)

    f = f1 * f2
    assert f.x == f1.x * f2.x
    assert f.func(3) == f1.func(3) * f2.func(3)

    f = max(f1, f2)
    assert f.x == max(f1.x, f2.x)
    assert f.func(3) == max(f1.func(3), f2.func(3))
