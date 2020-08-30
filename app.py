import streamlit as st
import SessionState as stss
import zipfile
import tempfile
import pprint
import os
import numpy as np
import matplotlib.pyplot as plt
from examples.HanLab_calculate_ODNP_FilterCorruptData import hanlab_calculate_odnp
from dnpLab.dnpHydration import HydrationParameter
from app_helper import ProcParameter, dict_to_str, get_table_download_link

st.set_option('deprecation.showfileUploaderEncoding', False)

print = pprint.pprint

# TEMPDIR = '/tmp/odnplab/'
TEMPDIR = None
VERSION = "v1.3"
CNSI_EMX_LINK = 'https://www.mrl.ucsb.edu/spectroscopy-facility/instruments/7-bruker-emxplus-epr-spectrometer'
DEMO_DATA_LINK = 'https://github.com/ylin00/odnplab/raw/master/20190821_TW_4OH-TEMPO_500uM_.zip'
ISSUE_COMPLAINT_LINK = 'https://github.com/ylin00/odnplab/issues'
DNPLAB_REPO_LINK = 'https://github.com/DNPLab/dnpLab'
DNPLAB_DOC_LINK = 'http://dnplab.net/'


def set_par(ppar:ProcParameter, hpar:HydrationParameter):
    """Prompt for users to choose parameters

    Returns: tuple(ProcParameter, HydrationParameter)

    """
    # Defaults
    ppar.eiw = 20
    hpar.field = 348.5
    hpar.t1_interp_method = 'second_order'

    # st.sidebar.markdown('**Experiments**')
    hpar.spin_C = st.sidebar.number_input(
        "Spin label concentration (uM)", min_value = 0.01, value=500.0, step=1.0, key='spin_C'
    )
    hpar.T100 = st.sidebar.number_input(
        "T1,0(0) (s)", min_value=0.1, max_value=3.0, value=2.5, step=0.05, key='t100'
    )
    hpar.smax_model = st.sidebar.radio(
        'The spin is ', options=['tethered', 'free'], key='smax_model'
    )

    if st.sidebar.checkbox("More"):
        ppar.eiw = st.sidebar.number_input(
            "Integration Width", min_value=10, max_value=500, value=20, step=10, key='eiw'
        )
        hpar.field = st.sidebar.number_input(
            "Field (mT)", value=348.5, step=1.0, key='field'
        )
        hpar.t1_interp_method = st.sidebar.radio(
            'T1 interpolation method', options=['linear', 'second_order'], index=1, key='t1_interp_method'
        )

    return ppar, hpar


def run(uploaded_file, ppar:ProcParameter, hpar:HydrationParameter):
    """Process uploaded zipfile

    Args:
        uploaded_file: zip file object

    Returns: tuple(dict, str, dict)
        mydict: dictionary of hydration results in strings
        expname: name of the experiment
        hresults: dictionary of hydration results

    """
    # print(f"You just upload this file -> {uploaded_file}")
    # print(f"But I am in a demo mode and not going to run it actually")

    with tempfile.TemporaryDirectory(dir=TEMPDIR) as tmpdir:

        # upzip to tmpdir
        with zipfile.ZipFile(uploaded_file, "r") as zip_ref:
            zip_ref.extractall(tmpdir)
            # Select the first folder ended with '/1/', no matter how deep
            expname = sorted([x for x in zip_ref.namelist() if x[-3:] == '/1/' and 'pdata' not in x])
            if expname is None or len(expname) == 0:
                st.warning(f"""
                I could not find a folder with experiment number 1. \n
                Could you double check if you have `my_odnp_exp/1/`? \n
                If problems are still there, please report the issue below.
                 """)
                return {}, '', {}
            else:
                expname = expname[0][0:-2]

        # Process CNSI ODNP and return a str of results
        path = os.path.join(tmpdir, expname)  # path to CNSI data folder
        pars = {
            'integration_width'  : ppar.eiw,
             'spin_C'             : hpar.spin_C,
             'field'              : hpar.field,
             'T100'               : hpar.T100,
             'smax_model'         : hpar.smax_model,
             't1_interp_method'   : hpar.t1_interp_method,
             'drop_e_powers'       : ppar['drop_e_powers'],
             'drop_t1_powers'      : ppar['drop_t1_powers']
        }  # TODO: creating a dictionary is error-prone, replace it with a parameter class
        hresults = hanlab_calculate_odnp(path, pars, verbose=ppar.verbose)
        # Check T1,0 vs T1,0,0
        t10, t10std, t100 = hresults['T10'], hresults['T10_std'], hpar.T100
        if t10 + t10std > t100:
            st.warning(
                r"Error: $T_{1,0,0}$ must no less than T_{1,0} + stdev(T_{1,0}) \n"+
                r"$T_{1,0,0}, T_{1,0}, stdev(T_{1,0}) = "+
                rf"{round(t100,2)}, {round(t10,2)}, {round(t10std,2)}$")
            return {}, '', {}
        mydict = {k: v for k, v in hresults.items()
                  if type(v) != type(np.ndarray([]))}
        mydict.update({k: ', '.join([f"{vi:.4f}" for vi in v])
                       for k, v in hresults.items()
                       if type(v) == type(np.ndarray([]))})

    return mydict, expname, hresults


def plot(data:dict):
    """Create a plot

        +-----------+
        |     A     |
        |           |
        +-----+-----+
        |  B  |  C  |
        +-----+-----+

        A: ksigma ~ power
        B: E ~ power
        C: T1 ~ power

    Args:
        data: a dictionary of hydration results

    """
    if not data:
        return

    def plot_t1(ax, x, y, label=None):
        ax.plot(x, y, color = '#003660', marker = 'o', linestyle = 'none', label=label)

    def plot_t1_fit(ax, x, y, label=None):
        ax.plot(x, y, color = '#F37021', label=label)

    def plot_enhancement(ax, x, y, label=None):
        ax.plot(x, y, color='#003660', marker='o', linestyle='none', label=label)

    def plot_enhancement_fit(ax, x, y, label='Fit'):
        plot_t1_fit(ax, x, y, label)

    def plot_ksigma(ax, x, y, label=r'DNPLab $k_\sigma$[p]'):
        plot_t1(ax, x, y, label)

    def plot_ksigma_fit(ax, x, y, label=r'dnpHydration Fit'):
        plot_t1_fit(ax, x, y, label)

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    f3_ax1 = fig.add_subplot(gs[:, 0:2])
    f3_ax2 = fig.add_subplot(gs[0, 2])
    f3_ax3 = fig.add_subplot(gs[1, 2])
    # ksigma plot
    plot_ksigma(f3_ax1, data['E_power'], data['ksigma_array'], label=None)
    plot_ksigma_fit(f3_ax1, data['E_power'], data['ksigma_fit'], label='Fit')
    f3_ax1.set_xlabel('Power')
    f3_ax1.set_ylabel(r'$k_\sigma$ ($s^{-1} M^{-1}$)')
    f3_ax1.legend()
    # text in ksigma plot
    x_max, y_max = max(data['E_power']), max(data['ksigma_array'])
    for offset, text in zip([0, 0.05, 0.10, 0.15, 0.2, 0.25], [
        rf"$k_\sigma = {round(data['ksigma'], 2)} \pm {round(data['ksigma_stdd'])}$" + r" $s^{-1} M^{-1}$",
        rf"$k_\rho = {round(data['krho'], 2)}$" + r" $s^{-1} M^{-1}$",
        r"$k_{low}"+rf"= {round(data['klow'], 2)}$" + r" $s^{-1} M^{-1}$",
        r"$t_{corr}"+rf"= {round(data['tcorr'], 2)}$ ps",
        r"$D_{local}"+rf"= {round(data['Dlocal']*1e9, 3)}$"+r"$\times10^{-9}d^2/s$",
        rf"$\xi = {round(data['coupling_factor'], 4)}$"
    ]):
        f3_ax1.text(x_max * 0.4, y_max * (0.5 - offset), text, fontsize=12)
    # Enhancement plot
    plot_enhancement(f3_ax2, data['E_power'], data['E'], label=None)
    plot_enhancement_fit(f3_ax2, data['E_power'], data['uncorrected_Ep'], label='Fit')
    f3_ax2.set_xlabel('Power')
    f3_ax2.set_ylabel('Enhancement')
    f3_ax2.legend()
    # T1 plot
    plot_t1(f3_ax3, data['T1_power'], data['T1'], label=None)
    plot_t1_fit(f3_ax3, data['E_power'], data['interpolated_T1'], label='Fit')
    f3_ax3.set_xlabel('Power')
    f3_ax3.set_ylabel(r'$T_1$ (s)')
    f3_ax3.legend()
    # text in T1 plot
    x_max = max(data['E_power'])
    y_min, y_max = min(data['interpolated_T1']), max(data['interpolated_T1'])
    f3_ax3.text(
        x_max * 0.15, y_min + (y_max-y_min) * 0.05,
        r"$T_{1,0}"+rf" = {round(data['T10'], 2)} \pm {round(data['T10_std'], 2)} s$"
    )
    st.pyplot(fig=fig)


def drop_data(drop_e_powers:list, drop_t1_powers:list):
    """Create selectbox for dropping bad data points

    Args:
        drop_t1_powers: list of T1 powers to choose from
        drop_e_powers: list of Enhancement powers to choose from

    Returns:
        Tuple(list, list): Tuple of selected list of enhancement powers and list of T1 powers

    """
    drop_es, drop_t1s = {}, {}
    if len(drop_e_powers) + len(drop_t1_powers) > 0:
        drop_es = st.sidebar.multiselect(
            'Drop Enhancements at power(s):', drop_e_powers, key='drop_es'
        )
        drop_t1s = st.sidebar.multiselect(
            'Drop T1 at power(s):', drop_t1_powers, key='drop_t1s'
        )
    return drop_es, drop_t1s


# =======THE APP=======
ss = stss.get(
    ppar = ProcParameter(drop_e_powers=[], drop_t1_powers=[]),
    hpar = HydrationParameter(),
    results = {},
    expname = '',
    old_expname='',
    data = {},
    epowers = [],   # E_power from the data
    t1powers = [],  # T1_power from the data
    b_run=False
)

st.markdown('<style>' + open('app.css').read() + '</style>', unsafe_allow_html=True)
st.title(f'ODNPLab: One-Step ODNP Processing \n {VERSION} \t Powered by [DNPLab]({DNPLAB_DOC_LINK}) ')

st.markdown("## Upload a Zip file")
uploaded_file = st.file_uploader("Here ->", type="zip")

if uploaded_file is not None:

    # Parameters
    ss.ppar.verbose = False
    ss.ppar, ss.hpar = set_par(ss.ppar, ss.hpar)

    # Process the data
    b_run = st.button("Run")
    if b_run:
        with st.spinner('This should take 10 seconds ...'):
            ss.results, ss.expname, ss.data = run(uploaded_file, ppar=ss.ppar, hpar=ss.hpar)
    else:
        st.markdown("^ Click Me ")

    # Present the results
    if b_run or ss.b_run:
        ss.b_run = True
        plot(ss.data)
        st.markdown(
            get_table_download_link(dict_to_str(ss.results), filename=ss.expname),
            unsafe_allow_html=True
        )
        st.write(ss.results)

        # Filter bad data points when results are present
        if ss.old_expname != ss.expname:
            ss.old_expname = ss.expname
            ss.epowers = ss.data['E_power']
            ss.t1powers = ss.data['T1_power']
        ss.ppar['drop_e_powers'], ss.ppar['drop_t1_powers'] = drop_data(ss.epowers, ss.t1powers)

st.markdown(f"""
## How to use
1. Collect your ODNP data on [UCSB CNSI EMXplus]({CNSI_EMX_LINK}).
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
- For windows 7 and above you can use 'add to zip file'.
- For Mac you can use 'compress'.

5. Upload the zip file and click run.
""")

st.markdown(f"""
## Demo
6. For demo, click [here]({DEMO_DATA_LINK}) to download a zip file and upload. The demo data came from (500 $\mu$M 4OH-TEMPO in water, {'$k_{sigma} = 95 s^{-1} M^{-1}$'}).

## Issues/Support
7. Report any issue [here]({ISSUE_COMPLAINT_LINK}) and I will get back to you shortly.
""")
