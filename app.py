import streamlit as st


def run(uploaded_file):
    print(f"You just upload this file -> {uploaded_file}")
    print(f"But I am in a demo mode and not going to run it actually")
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

