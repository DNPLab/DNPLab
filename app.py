import streamlit as st
import numpy as np
import zipfile
import tempfile
import pprint
import os
import csv
import base64
from odnpLab.examples.workupCNSI import ProcParameter, process_cnsi
from odnpLab.hydration import HydrationParameter, HydrationCalculator

print = pprint.pprint

# TEMPDIR = '/tmp/odnplab/'
TEMPDIR = None


def get_table_download_link(temp_file_path, filename='results'):
    """Generates a link allowing a temp_file_path to be downloaded
    in:  dataframe
    out: href string
    """
    b64 = base64.b64encode(temp_file_path.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}.txt">Download Results</a>'
    return href


def set_ppar(ppar:ProcParameter):
    """Prompt for users to choose parameters

    Returns: tuple(ProcParameter, HydrationParameter)

    """
    # Processing parameter
    st.sidebar.markdown("**Enhancement integration**")
    ppar.eic = st.sidebar.slider("Center", min_value=-100, max_value=100, value=0, step=10, key='eic')
    ppar.eiw = st.sidebar.slider("Width", min_value=10, max_value=500, value=250, step=10, key='eiw')
    st.sidebar.markdown("**T1 integration**")
    if not st.sidebar.checkbox("Same as Enhancement", value=True):
        ppar.tic = st.sidebar.slider("Center", min_value=-100, max_value=100, value=0, step=10, key='tic')
        ppar.tiw = st.sidebar.slider("Width", min_value=10, max_value=500, value=250, step=10, key='tiw')
    else:
        ppar.tic = ppar.eic
        ppar.tiw = ppar.eiw

    return ppar


def set_hpar(hpar:HydrationParameter):
    """Prompt for users to choose parameters

    Returns: tuple(ProcParameter, HydrationParameter)

    """
    # Processing parameter
    st.sidebar.markdown('**Hydration**')
    hpar.spin_C = st.sidebar.number_input("Spin label concentration (uM)", value=500.0, step=1.0, key='spin_C')
    hpar.field = st.sidebar.number_input("Field (mT)", value=hpar.field, step=1.0, key='field')
    hpar.smax_model = st.sidebar.radio('The spin is ', options=['tethered', 'free'], key='smax_model')
    hpar.t1_interp_method = st.sidebar.radio('T1 interpolation method', options=['linear', 'second_order'], index=0, key='t1_interp_method')

    # if st.sidebar.button('More'):
    #     st.write("No more")

    return hpar


def run(uploaded_file, ppar:ProcParameter, hpar:HydrationParameter):
    """

    Args:
        uploaded_file: zip file object

    Returns: tuple(dict, str)
        mydict: dictionary of results
        expname: name of the experiment

    """
    # print(f"You just upload this file -> {uploaded_file}")
    # print(f"But I am in a demo mode and not going to run it actually")

    with tempfile.TemporaryDirectory(dir=TEMPDIR) as tmpdir:

        # upzip to tmpdir
        with zipfile.ZipFile(uploaded_file, "r") as zip_ref:
            zip_ref.extractall(tmpdir)
            expname = zip_ref.namelist()[0]

        # Process CNSI ODNP and return a str of results
        path = os.path.join(tmpdir, expname)  # path to CNSI data folder

        rest = process_cnsi(path, ppar)

        t1, t1_power, e, e_power = rest['T1p'], rest['T1powers'], rest['Ep'], rest['Epowers']

        hpar.T10 = rest['T10']

        hc = HydrationCalculator(T1=t1, T1_power=t1_power, E=e, E_power=e_power, hp=hpar)
        hc.run()

        # Create dictionary of results
        mydict = {k:v for k, v in hc.results.__dict__.items()
                  if type(v) != type(np.ndarray([]))}
        mydict.update({k: ', '.join([f"{vi:.4f}" for vi in v])
                       for k, v in hc.results.__dict__.items()
                       if type(v) == type(np.ndarray([]))})

    return mydict, expname


def dict_to_str(mydict):
    mylist = [f"{k} \t {v}" for k, v in mydict.items()]
    return '\n'.join(mylist)


# =======THE APP=======
st.title('ODNPLab \n One-step ODNP processing')
st.markdown("""
## How to use

1. Collect your ODNP data in UCSB CNSI (talk to the facility)
2. Save your data in an experiment folder. For demo only here we use `my_odnp_exp`.
3. Your experiment folder should look like the following:
```
my_odnp_exp/
            1/...
            2/...
            3/...
            ...
            t1_powers.mat
            power.mat
```
4. Right click the experiment folder `my_odnp_exp` and create a zip file:
- For windows 7 and above you can use 'add to zip file'
- For Mac you can use 'compress'

5. Upload the zip file and click run
""")

# 6. Download a demo (500 $\mu$M 4OH-TEMPO in water, $k_{sigma} ~= 95 s^{-1} M^{-1}$)
# st.write(f'<a href="" download="ucsb_cnsi_odnp_demo.zip">Download Results</a>')
# _ = st.text_input('Your Lucky Number', value='42')

st.markdown("## Upload a Zip file")
uploaded_file = st.file_uploader("Here ->", type="zip")

if uploaded_file is not None:

    # Parameters
    ppar = ProcParameter()
    ppar.verbose = False
    ppar = set_ppar(ppar)

    hpar = HydrationParameter()
    hpar = set_hpar(hpar)

    if st.button("Run"):

        results, expname = run(uploaded_file, ppar=ppar, hpar=hpar)
        st.write(results)
        st.markdown(get_table_download_link(dict_to_str(results), filename=expname),
                    unsafe_allow_html=True)
    else:
        st.write("^ Click Me ")
