class FragmentRetroError(Exception):
    """Base class for all exceptions in FragmentRetro."""

    pass


class SmartsParsingError(FragmentRetroError, ValueError):
    """Error raised when a SMARTS string is parsed incorrectly."""

    pass
