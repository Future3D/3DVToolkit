from inspect import signature, Parameter


def crosscheck_kwargs_function(kwargs, function):
    f_positional = []
    f_kwargs = []
    # accepts_args = False
    accept_kwargs = False
    for x, p in signature(function).parameters.items():
        if p.default == Parameter.empty:
            if p.kind == Parameter.VAR_POSITIONAL:
                # accepts_args = True
                pass
            elif p.kind == Parameter.VAR_KEYWORD:
                accept_kwargs = True
            else:
                f_positional.append(x)
        else:
            f_kwargs.append(x)

    valid_kwargs = {}

    for x in f_positional:
        valid_kwargs[x] = kwargs[x]
        del kwargs[x]

    for x in f_kwargs:
        if x in kwargs:
            valid_kwargs[x] = kwargs[x]
            del kwargs[x]

    if accept_kwargs:
        for x in kwargs:
            valid_kwargs[x] = kwargs[x]

    return valid_kwargs
