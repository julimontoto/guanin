![Logo](/image/logoguanin.png)
  
# GUANIN

| GUi-driven Analyser for Nanostring Interactive Normalization |

[Explore the docs](https://github.com/julimontoto/guanin)

[View Demo](https://github.com/github_username/repo_name) TBD

[Report But](https://github.com/julimontoto/guanin/issues)

[Request Feature](https://github.com/github_username/repo_name/issues)

## Table of Contents

* [About Guanin](#about-guanin)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgments](#acknowledgements)


![Guanin Screen Shot](https://i.imgur.com/TBTcTnm.png)

## Getting started

### Prerequisites

GUANIN should run under Linux, MacOS and Windows. You'll need Python >= 3.9 installed.

### Installation

Check INSTALL.txt for further details.

## Running

### The CLI

Open a console and from any path, run the following:

    $ guanin-cli

### The GUI

A simple GUI is included using [pyQT6](https://pypi.org/project/PyQt6/).

From a console, launch the GUI with:

python3 guanin-gui.py

This *executable* is installed in your PATH (TBD!). You can find it and launch
it with a double-click.
Although, running it through the console will provide further information about the execution of the program.


<!-- LICENSE -->
## License

Distributed under the GPL License. See `LICENSE.txt` for more information.

<!-- For developers -->
## For developers

Pull requests are only welcome if both the following are true:

1. They include at least some test.

2. They solve a bug.

We don't accept pull requests of enhancements, new requirements, "It would be
nice if...", and so on. If that is your case, fork and develop.

### For internal usage only

To create a new release, follow the steps:

1. Bump the version at `nqcview/__init__.py`.???

2. Create a new bundle **locally** with `python -m build`.

3. Check you can install the bundle **locally**, in a fresh Virtualenv, with:

    `$ pip install NanostringQC-0.1.0.tar.gz`????

    **AND**

     `$ pip install NanostringQC-0.1.0-py3-none-any.whl`????

4. Check it's working as expected, both the CLI and the GUI.

5. Bump the version in this Readme, in the install link above.

6. Put the code under the CVS, **tag it with version**, and push.

7. At Github, create a new Release, selecting the tag from above, and attach
  the binaries tested in [3].


## Contact

Julián Montoto-Louzao - [@julimontoto](https://twitter.com/julimontoto) - juli.mlouzao@gmail.com

GUANIN: [https://github.com/julimontoto/guanin](https://github.com/julimontoto/guanin)
