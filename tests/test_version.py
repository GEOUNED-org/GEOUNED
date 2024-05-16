import geouned


def test_version():
    """Check that __version__ exists and is correctly formatted"""
    version = geouned.__version__
    # Ensure it is given as a string
    assert isinstance(version, str)

    # develop is used as the default version name in the dockerfile
    if version != "develop":
        version_bits = version.split(".")
        assert len(version_bits) >= 2, f"Version number is {version}"
        # Ensure both of the first two parts is convertable to int
        # (The third/fourth part, if it exists, may be a development version)
        for bit in version_bits[:2]:
            assert bit.isdigit(), f"Version number is {version}"
