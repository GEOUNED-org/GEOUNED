import string

version = "3.6"


def make_label(item):
    if 0 < len(item):
        return "{" + item + "}"
    else:
        return item


class SafeDict(dict):

    def __getitem__(self, item):
        return super(SafeDict, self).__getitem__(item) or ""

    def __missing__(self, key):
        return make_label(key)


# ' 1 #'
class PartialFormatter(string.Formatter):

    def format(*args, **kwargs):
        if not args:
            raise TypeError(
                "descriptor 'format' of 'Formatter' object " "needs an argument"
            )
        self, args = args[0], args[1:]  # allow the "self" keyword be passed
        try:
            format_string, args = (
                args[0],
                args[1:],
            )  # allow the "format_string" keyword be passed
        except IndexError:
            if "format_string" in kwargs:
                format_string = kwargs.pop("format_string")
            else:
                raise TypeError(
                    "format() missing 1 required positional "
                    "argument: 'format_string'"
                )
        return string.Formatter.vformat(self, format_string, args, SafeDict(kwargs))

    def _vformat(
        self, format_string, args, kwargs, used_args, recursion_depth, auto_arg_index=0
    ):
        """Clone from Python3 version to fix Python2 mess with unnamed {} format specifiers"""
        if recursion_depth < 0:
            raise ValueError("Max string recursion exceeded")
        result = []
        for literal_text, field_name, format_spec, conversion in self.parse(
            format_string
        ):

            # output the literal text
            if literal_text:
                result.append(literal_text)

            # if there's a field, output it
            if field_name is not None:
                # this is some markup, find the object and do
                #  the formatting

                # handle arg indexing when empty field_names are given.
                if field_name == "":
                    if auto_arg_index is False:
                        raise ValueError(
                            "cannot switch from manual field "
                            "specification to automatic field "
                            "numbering"
                        )
                    field_name = str(auto_arg_index)
                    auto_arg_index += 1
                elif field_name.isdigit():
                    if auto_arg_index:
                        raise ValueError(
                            "cannot switch from manual field "
                            "specification to automatic field "
                            "numbering"
                        )
                    # disable auto arg incrementing, if it gets
                    # used later on, then an exception will be raised
                    auto_arg_index = False

                # given the field_name, find the object it references
                #  and the argument it came from
                obj, arg_used = self.get_field(field_name, args, kwargs)
                used_args.add(arg_used)

                # do any conversion on the resulting object
                obj = self.convert_field(obj, conversion)

                # expand the format spec, if needed
                if version == "2.7":
                    format_spec = self._vformat(
                        format_spec,
                        args,
                        kwargs,
                        used_args,
                        recursion_depth - 1,
                        auto_arg_index=auto_arg_index,
                    )
                else:
                    format_spec, auto_arg_index = self._vformat(
                        format_spec,
                        args,
                        kwargs,
                        used_args,
                        recursion_depth - 1,
                        auto_arg_index=auto_arg_index,
                    )
                # format the object and append to the result
                result.append(self.format_field(obj, format_spec))

        result = "".join(result)

        if version == "2.7":
            return result
        else:
            return result, auto_arg_index
