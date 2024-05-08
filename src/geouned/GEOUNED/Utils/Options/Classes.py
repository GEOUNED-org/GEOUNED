class Options:
    from .optionsDefault import default_values, type_dict

    @classmethod
    def set_default_attribute(cls):
        for key, value in cls.default_values.items():
            setattr(cls, key, value)

    @classmethod
    def set_attribute(cls, key, value):
        if key in cls.default_values.keys():
            setattr(cls, key, value)


class Tolerances:
    from .tolerancesDefault import KwrdEquiv, default_values, type_dict

    @classmethod
    def set_default_attribute(cls):
        for key, value in cls.default_values.items():
            setattr(cls, key, value)

    @classmethod
    def set_attribute(cls, key, value):
        if key in cls.default_values.keys():
            setattr(cls, key, value)


class McnpNumericFormat:
    from .mcnpNumericDefault import default_values

    @classmethod
    def set_default_attribute(cls):
        for key, value in cls.default_values.items():
            setattr(cls, key, value)

    @classmethod
    def set_attribute(cls, key, value):
        if key in cls.default_values.keys():
            setattr(cls, key, value)
