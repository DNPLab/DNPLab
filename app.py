import streamlit as st
import zipfile
import tempfile
import pprint
import os
from dnpLab.examples.workupCNSI import ProcParameter,
print = pprint.pprint

# TEMPDIR = '/tmp/odnplab/'
TEMPDIR = None

def run(uploaded_file):
    print(f"You just upload this file -> {uploaded_file}")
    print(f"But I am in a demo mode and not going to run it actually")

    # upzip
    with zipfile.ZipFile(uploaded_file, "r") as zip_ref:
        # print(f"zipfile filename {zip_ref.filename}")
        # print(f"zipfile filelist {zip_ref.filelist}")
        # print(f"zipfile namelist {zip_ref.namelist()}")
        odnpdir = tempfile.mkdtemp(dir = TEMPDIR)

        # After extraction, a lot of OS unrelated folders will be there
        zip_ref.extractall(odnpdir)

        odnppath = os.path.join(odnpdir, zip_ref.namelist()[0])
        print(f"odnppath= {odnppath}")


    # process

    # return

    pass

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
        run(uploaded_file)
    else:
        st.write("^ Click Me ")

