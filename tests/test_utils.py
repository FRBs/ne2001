from ne2001.utils import ClassOperation


class Foo(ClassOperation):
    """
    """

    def __init__(self, func, x):
        self.x = x
        self.func = func

def test_classoperation():
    f1 = Foo(lambda x: 2*x, 1)
    f2 = Foo(lambda x: x**2, 2)
    f3 = Foo(lambda x: x**3, 3)

    f = f1 + f2 + f3
    assert f.x == f1.x + f2.x + f3.x
    assert f.func(3) == f1.func(3) + f2.func(3) + f3.func(3)

    f = f1 - f2 - f3
    assert f.x == f1.x - f2.x - f3.x
    assert f.func(3) == f1.func(3) - f2.func(3) - f3.func(3)

    f = f1 * f2 * f3
    assert f.x == f1.x * f2.x * f3.x
    assert f.func(3) == f1.func(3) * f2.func(3) * f3.func(3)

    f = max(max(f1, f2), f3)
    assert f.x == max(max(f1.x, f2.x), f3.x)
    assert f.func(3) == max(max(f1.func(3), f2.func(3)), f3.func(3))
