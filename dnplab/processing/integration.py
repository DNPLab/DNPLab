def integrate(
    all_data,
    dim="f2",
    type="trapz",
    integrate_center=0,
    integrate_width="full",
):
    """Integrate data down given dimension

    Args:
        all_data (dnpdata,dict): Data container

    +------------------+---------------+----------+-------------------------------+
    | parameter        | type          | default  | description                   |
    +==================+===============+==========+===============================+
    | dim              | str           | 'f2'     | dimension to integrate        |
    +------------------+---------------+----------+-------------------------------+
    | type             | str           | 'single' | 'single' or 'double' integral |
    +------------------+---------------+----------+-------------------------------+
    | integrate_center | float or list | 0        | center of integration window  |
    +------------------+---------------+----------+-------------------------------+
    | integrate_width  | float or list | "full"   | width of integration window   |
    +------------------+---------------+----------+-------------------------------+

    Returns:
        dnpdata: integrals of data
    """

    data, isDict = return_data(all_data)
    index = data.dims.index(dim)

    data_new = data.copy()
    if type == "double":
        first_int = scipy.integrate.cumtrapz(
            data.values, x=data.coords[dim], axis=index, initial=0
        )
        data_new.values = first_int

    if integrate_width == "full":
        pass
    elif isinstance(integrate_width, (int, float)) and isinstance(
        integrate_center, (int, float)
    ):
        integrateMin = integrate_center - np.abs(integrate_width) / 2.0
        integrateMax = integrate_center + np.abs(integrate_width) / 2.0
        data_new = data_new[dim, (integrateMin, integrateMax)]

    elif (
        isinstance(integrate_width, list)
        and isinstance(integrate_center, list)
        and all((isinstance(x, (int, float)) for x in integrate_width))
        and all((isinstance(x, (int, float)) for x in integrate_center))
    ):
        if len(integrate_width) != len(integrate_center):
            raise TypeError(
                "If integrate_center and integrate_width are both lists, they must be the same length"
            )

        integrateMinMax = [
            [
                (cent - np.abs(integrate_width[x]) / 2.0),
                (cent + np.abs(integrate_width[x]) / 2.0),
            ]
            for x, cent in enumerate(integrate_center)
        ]
        data_new = [data_new[dim, (mx[0], mx[1])] for mx in integrateMinMax]
    elif (
        isinstance(integrate_width, (int, float))
        and isinstance(integrate_center, list)
        and all((isinstance(x, (int, float)) for x in integrate_center))
    ):
        integrateMinMax = [
            [
                (cent - np.abs(integrate_width) / 2.0),
                (cent + np.abs(integrate_width) / 2.0),
            ]
            for cent in integrate_center
        ]
        data_new = [data_new[dim, (mx[0], mx[1])] for mx in integrateMinMax]

    elif (
        isinstance(integrate_center, (int, float))
        and isinstance(integrate_width, list)
        and all((isinstance(x, (int, float)) for x in integrate_width))
    ):
        integrateMinMax = [
            [
                (integrate_center - np.abs(wid) / 2.0),
                (integrate_center + np.abs(wid) / 2.0),
            ]
            for wid in integrate_width
        ]
        data_new = [data_new[dim, (mx[0], mx[1])] for mx in integrateMinMax]

    else:
        raise ValueError(
            "integrate_width must be 'full', int, float, or list of int or float; integrate_center must be int, float, or list of ints or floats"
        )

    remaining_dims = [x for x in data.dims if x != dim]
    if (
        len(remaining_dims) == 0
        and isinstance(integrate_center, (int, float))
        and isinstance(integrate_width, (int, float))
    ):
        remaining_dims = ["index"]
        remaining_coords = [np.array([0])]
    else:
        remaining_coords = [data.coords[x] for x in remaining_dims]

    if isinstance(data_new, list):
        data_integrals = [
            np.trapz(x.values, x=x.coords[dim], axis=index) for x in data_new
        ]
        if not all([isinstance(x, list) for x in [integrate_center, integrate_width]]):
            if isinstance(integrate_center, list) and not isinstance(
                integrate_width, list
            ):
                remaining_coords = [np.array(integrate_center)] + remaining_coords
                remaining_dims = ["center"] + remaining_dims
            elif isinstance(integrate_width, list) and not isinstance(
                integrate_center, list
            ):
                remaining_coords = [np.array(integrate_width)] + remaining_coords
                remaining_dims = ["width"] + remaining_dims
            data_values = np.array(data_integrals)
        elif isinstance(integrate_center, list) and isinstance(integrate_width, list):
            ind_dim = list(set(data.dims) - set([dim]))[0]
            ind_shape = data.shape[data.index(ind_dim)]
            remaining_coords = [
                np.array(integrate_center),
                np.array(integrate_width),
            ] + remaining_coords
            remaining_dims = ["center", "width"] + remaining_dims
            data_values = np.array(
                tuple([data_integrals for _ in range(len(data_integrals))])
            ).reshape(len(data_integrals), len(data_integrals), ind_shape)
    else:
        data_values = np.trapz(data_new.values, x=data_new.coords[dim], axis=index)

    if not isinstance(data_values, (list, np.ndarray)):
        data_values = [data_values]

    integrate_data = dnpdata(np.array(data_values), remaining_coords, remaining_dims)

    integrate_data.attrs["integrate_center"] = integrate_center
    integrate_data.attrs["integrate_width"] = integrate_width

    if type == "double":
        integrate_data.attrs["first_integral"] = first_int
        integrate_data.attrs["dim_coords"] = data.coords[dim]

    if isDict:
        all_data["integrals"] = integrate_data
    else:
        return integrate_data