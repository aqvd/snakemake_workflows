def field_from_sample(Sample, field):
    # ?P<group_name> for clarity
    match = re.match(
        r"(?P<readID>\w+)_(?P<readLetter>\w+)_(?P<readRun>\w+)_(?P<readLane>\w+)", Sample)

    return match.group(field)

def get_column_df(df , column, filt_out, **kwags):
    """
    Get column, eg TumorSample or NormalSample, for individual.
    Optional filter value for the results
    """
    # Filter using **kwags
    row_filters = []
    for k,v in kwags.items():
        row_filters.append('({0} == "{1}")'.format(k, v))

    query_expr = " & ".join(row_filters)

    # query expression and Select column with .loc[:, column]
    res = df.query(query_expr).loc[:, column]

    # Filter values of that column
    res = np.array(list(filter(lambda x: x != filt_out, res)))
    return(res)

def repeat_argument(argument, value_list, **kwags):
    """
    sometimes a argument (eg: -I <file>) must be specified multiple times,
    this function creates the appropriate concatenation.
    Include whitespace in argument (eg "-I ") if needed
    """

    argument_l = [str(argument) + str(v) + " " for v in value_list]
    res = "".join(argument_l)

    return(res)

def expand_argument(path, string_to_expand,
                    df, column, filt_out,
                    argument, **kwags):
    """
    generates repeated argument value pairs (eg -I path/file_1.txt -I path/file_2.txt)
    using expand function from snakemake to create the paths.

    string_to_expand must be '{expansion}_bla_bla_bla.xxx'
    """
    if not path.endswith("/"):
        path = path + "/"

    # get the values from df
    to_expand = get_column_df(df, column, filt_out, **kwags)

    # expand them to file paths
    arg_value_l = expand(path + string_to_expand, expansion= to_expand)

    # Join using the argument
    res = repeat_argument(argument, arg_value_l)

    return(res)