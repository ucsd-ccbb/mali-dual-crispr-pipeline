def split_delimited_string_to_list(input_str, cast_func=str, delimiter_str=","):
    return [cast_func(x.strip()) for x in input_str.split(delimiter_str)]