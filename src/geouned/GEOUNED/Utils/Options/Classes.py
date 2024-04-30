class Options:
    from .optionsDefault import defaultValues, typeDict

    @classmethod
    def setDefaultAttribute(cls):
        for key, value in cls.defaultValues.items():
            setattr(cls, key, value)

    @classmethod
    def setAttribute(cls, key, value):
        if key in cls.defaultValues.keys():
            setattr(cls, key, value)


class Tolerances:
    from .tolerancesDefault import defaultValues, typeDict, KwrdEquiv

    @classmethod
    def setDefaultAttribute(cls):
        for key, value in cls.defaultValues.items():
            setattr(cls, key, value)

    @classmethod
    def setAttribute(cls, key, value):
        if key in cls.defaultValues.keys():
            setattr(cls, key, value)


class McnpNumericFormat:
    from .mcnpNumericDefault import defaultValues

    @classmethod
    def setDefaultAttribute(cls):
        for key, value in cls.defaultValues.items():
            setattr(cls, key, value)

    @classmethod
    def setAttribute(cls, key, value):
        if key in cls.defaultValues.keys():
            setattr(cls, key, value)
