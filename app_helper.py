import base64

from dnpLab.dnpHydration import Parameter


class ProcParameter(Parameter):
    """Processing Parameters

    Attributes:
        eic (float): enhancement data integration window center
        eiw (float): enhancement data integration window width
        tic (float): T1 data integration window center
        tiw (float): T1 data integration window width
        verbose (bool): Whether verbose.

    """
    def __init__(self, *args, **kwargs):  # TODO: enable manual adjustment
        super().__init__(*args, **kwargs)
        eic = 0
        eiw = 100
        tic = 0
        tiw = 100
        self.eic, self.eiw, self.tic, self.tiw = eic, eiw, tic, tiw
        self.verbose = True


def dict_to_str(mydict):
    mylist = [f"{k} \t {v}" for k, v in mydict.items()]
    return '\n'.join(mylist)


def get_table_download_link(temp_file_path, filename='results'):
    """Generates a link allowing a temp_file_path to be downloaded

    Args:
        temp_file_path(str): A string to write to a txt file and download.
        filename(str): the txt file name to generate.

    """
    b64 = base64.b64encode(temp_file_path.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}.txt">Download Results</a>'
    return href