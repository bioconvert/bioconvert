from bioconvert.core.decorators import requires


def test_require_binaries():
    def f():
        pass

    g = requires(external_binary="mv")(f)
    assert g.is_disabled is False
    # test cache now
    g = requires(external_binary="mv")(f)
    assert g.is_disabled is False
    g = requires(external_binaries=["ls", ])(f)
    assert g.is_disabled is False
    g = requires(external_binaries=["ls", "mv", ])(f)
    assert g.is_disabled is False

    g = requires(external_binary="tagada1")(f)
    assert g.is_disabled
    # test cache now
    g = requires(external_binary="tagada1")(f)
    assert g.is_disabled
    g = requires(external_binaries=["rm", "tagada2"])(f)
    assert g.is_disabled


def test_require_libraries():
    def f():
        pass

    g = requires(python_library="pandas")(f)
    assert g.is_disabled is False
    # test cache now
    g = requires(python_library="pandas")(f)
    assert g.is_disabled is False
    g = requires(python_libraries=["pip", ])(f)
    assert g.is_disabled is False
    g = requires(python_libraries=["pip", "pandas", ])(f)
    assert g.is_disabled is False

    g = requires(python_library="tagada3")(f)
    assert g.is_disabled
    # test cache now
    g = requires(python_library="tagada3")(f)
    assert g.is_disabled
    g = requires(python_libraries=["urllib3", "tagada4"])(f)
    assert g.is_disabled


def test_require_both():
    def f():
        pass

    g = requires(python_library="pandas", external_binary="top")(f)
    assert g.is_disabled is False

    g = requires(python_library="tagada5", external_binary="top")(f)
    assert g.is_disabled

    g = requires(python_library="pandas", external_binary="tagada6")(f)
    assert g.is_disabled

    g = requires(python_library="tagada7", external_binary="tagada8")(f)
    assert g.is_disabled