from bioconvert.core.decorators import requires


def test_require_binaries():
    def f():
        pass

    g = requires(external_binary="mv")(f)
    assert g.is_disabled is False
    g = requires(external_binaries=["ls", ])(f)
    assert g.is_disabled is False

    g = requires(external_binary="tralala")(f)
    assert g.is_disabled
    g = requires(external_binaries=["rm", "tagada"])(f)
    assert g.is_disabled


def test_require_libraries():
    def f():
        pass

    g = requires(python_library="csv")(f)
    assert g.is_disabled is False
    g = requires(python_libraries=["sys", ])(f)
    assert g.is_disabled is False

    g = requires(python_library="tralalaaa")(f)
    assert g.is_disabled
    g = requires(python_libraries=["os", "tagadazz"])(f)
    assert g.is_disabled


def test_require_both():
    def f():
        pass

    g = requires(python_library="csv", external_binary="top")(f)
    assert g.is_disabled is False

    g = requires(python_library="tralalaaazzz", external_binary="top")(f)
    assert g.is_disabled

    g = requires(python_library="csv", external_binary="topazeazeazeazeaze")(f)
    assert g.is_disabled

    g = requires(python_library="tralalaaazzzzz", external_binary="tagadaZZ")(f)
    assert g.is_disabled
