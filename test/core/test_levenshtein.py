from bioconvert.core.levenshtein import wf_levenshtein, classic_levenshtein


def test_wf_levenshtein():
    levenshtein_tests(classic_levenshtein)


def test_classic_levenshtein():
    levenshtein_tests(wf_levenshtein)


def levenshtein_tests(levenshtein):
    assert 1 == levenshtein("kitten", "kittenn")
    assert 3 == levenshtein("kitten", "sitting")
    assert 0 == levenshtein("kitten", "kitten")
    assert 0 == levenshtein("", "")
    assert 2 == levenshtein("sitting", "sititng")
