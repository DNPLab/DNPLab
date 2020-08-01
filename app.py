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


def run(uploaded_file):
    print(f"You just upload this file -> {uploaded_file}")
    print(f"But I am in a demo mode and not going to run it actually")

    # upzip
    with zipfile.ZipFile(uploaded_file, "r") as zip_ref:
        odnpdir = tempfile.mkdtemp(dir = TEMPDIR)
        # After extraction, a lot of OS unrelated folders will be there
        zip_ref.extractall(odnpdir)
        expname = zip_ref.namelist()[0]
        odnppath = os.path.join(odnpdir, expname)

    """Process CNSI ODNP and return a str of results"""
    path = odnppath  # path to CNSI data folder
    ppar = ProcParameter()
    hpar = HydrationParameter()

    ppar.verbose = False

    rest = process_cnsi(path, ppar)
    t1, t1_power, e, e_power = rest['T1p'], rest['T1powers'], rest['Ep'], rest['Epowers']
    print(t1_power)

    hpar.slC = 500
    hpar.T10 = rest['T10']
    hpar.field = 380.4

    hc = HydrationCalculator(T1=t1, T1_power=t1_power, E=e, E_power=e_power, hp=hpar)
    hc.run()
    print(hc.results)

    mydict = {k:v for k, v in hc.results.__dict__.items() if type(v) != type(np.ndarray([]))}
    mydict.update({k: ', '.join([f"{vi:.4f}" for vi in v]) for k, v in hc.results.__dict__.items() if
                   type(v) == type(np.ndarray([]))})

    return mydict, expname


def dict_to_str(mydict):
    mylist = [f"{k} \t {v}" for k, v in mydict.items()]
    return '\n'.join(mylist)

#
#
#
#
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
# _ = st.text_input('Your Lucky Number', value='42')

st.markdown("## Upload a Zip file")
uploaded_file = st.file_uploader("Here ->", type="zip")
if uploaded_file is not None:
    if st.button("Run"):
        results, expname = run(uploaded_file)
        st.write(results)
        st.markdown(get_table_download_link(dict_to_str(results), filename=expname),
                    unsafe_allow_html=True)
    else:
        st.write("^ Click Me ")

