# Contributing to DNPLab 
We love your input! We want to make contributing to this project as easy
and transparent as possible. The following contributions equally counts:

- Reporting a bug
- Proposing new features
- Improving documentation
- Adding a test case
- Submitting a fix
- Becoming a maintainer

## Bug reports and enhancement requests
We use GitHub issues to track public bugs and enhancement requests:
[opening a new issue](https://github.com/DNPLab/DNPLab/issues/new).

To report a bug, please:

1. Include a short, self-contained Python snippet reproducing the
   problem.
2. Explain why the current behavior is wrong/not desired and what you
   expect instead.


## Working with the code

### Git and Github

We will use git for version control, and github to host the DNPLab code
base.

Some resources for learning Git:
- [Atlassian Git tutorial](https://www.atlassian.com/git/tutorials/what-is-version-control)
- [Pandas's contribution guidance](https://pandas.pydata.org/pandas-docs/stable/development/contributing.html#version-control-git-and-github)

Some resources for learning GitHub:
- [Getting started with GitHub](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github)

### Forking and cloning

After [creating your free github account](https://github.com/join), you
can go to the [DNPLab project page](https://github.com/DNPLab/DNPLab)
and click the **Fork** button to creating a copy of the project to your
own Github account.

From now on, we will type commands in terminal if you are using MacOS or
Linux, or in cmd prompt if you are using windows (For Windows 10, make
sure you change to a directory with write permission).

After forking, you want to copy the project from Github to your local
machine by
```
git clone https://github.com/<your-github-user-name>/dnplab dnplab-yourname
cd dnplab-yourname
git remote add upstream https://github.com/dnplab/dnplab.git
```

The `master` branch only hosts the latest release, so you want to be on
the `develop` branch for the latest production-ready development. Type
```
git checkout develop
git pull upstream develop
```

### Creating a development environment

If you are making documentation changes only, you can skip this section.

To test out code changes, you need to set up all the dependencies.

We strongly recommend using virtual environment to avoid messing up with
your existing python environments. You can use conda by doing the
following
1. Install either
   [Anaconda](https://www.anaconda.com/products/individual) or
   [miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Make sure conda is up-to-date
   ```
   conda update conda
   ```
3. Create environment with Python>=3.6
   ```
   conda create --name environment_name python=3.8
   conda activate environment_name
   ```
4. Install dependencies. First, use `git status` to make sure you are at
   `develop` branch, then type
   ```
   python -m pip install -r requirements.txt
   ```

### Branching, commiting and pushing
When creating new branch, you want to make sure the `develop` is up to
date with DNPLab project. You can do so by
```
git checkout develop
git pull upstream develop
```

You want to branch out from the `develop` when making any change.
```
git branch yourname-gh-##
git checkout yourname-gh-##
```
We strongly recommend using `yourname-gh-##` as the branch name, where
`##` is the corresponding issue or pull request number in the DNPLab
project that you are addressing.

Now you are free to make any change to address corresponding issue. We
strongly recommend breaking tens of lines of changes into multiple
batches of small changes to ease the reviewing process.

After making changes, run the following and make sure no errors pop up
to ensure you are not breaking the code.
```
python -m pytest
```

Then run the following to check syntax error and correct formats
```
python -m flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
python -m black .
```

After all done. Commit your changes following
[pandas: committing your code](https://pandas.pydata.org/pandas-docs/stable/development/contributing.html#committing-your-code)

Push to your Github repository by
```
git push -u origin yourname-gh-##
```

### Pull requests
After reviewing your changes, you can file a pull request for the
maintainer of DNPLab to review and approve. See
[Github creating-a-pull-request](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request)

Make sure you are requesting to merge to `DNPLab/develop` from
`your-github-username/yourname-gh-##`.

Your changes will trigger multiple automatic checkings to ensure it
won't break the package. Then a maintainer from the DNPLab team will
accept or provide revising comments.

## Contributing to documentation


## License
By contributing, you agree that your contributions will be licensed under its MIT License.


## Becoming a maintainer

We strongly welcome and encourage committed individuals to help maintain
DNPLab. You don't have to be a python expert or a DNP expert. Please
contact any of the [current maintainers](http://dnplab.net/) if you are interested.
