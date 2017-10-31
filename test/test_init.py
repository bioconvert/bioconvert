from bioconvert.core.init import InitConverter



def test_init():
    init = InitConverter("bam", "sam")
    assert "BAM2SAM" in init.get_content()
