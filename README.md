
<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/julimontoto/guanin">
    <img src=![image/logoguanin.png]("image/logoguanin.png") width="156" height="156">
  </a>

<h3 align="center">GUANIN</h3>

  <p align="center">
    | GUi-driven Analyser for Nanostring Interactive Normalization |
    <br />
    <a href="https://github.com/julimontoto/guanin"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT GUANIN -->


[![Product Name Screen Shot][product-screenshot]]([https://example.com](https://i.imgur.com/TBTcTnm.png))

Aquí hay una captura representativa de GUANIN (que no se ve)

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->

### Prerequisites

GUANIN should run under Linux, MacOS and Windows.

### Installation

Assuming you have a Python >= 3.9 installed:
   
    $ pip install --user https://github.com/julimontoto/guanin/releases/download/0.1.0/GUANIN-0.1.0.tar.gz

Check INSTALL.txt for further details.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->

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

<p align="right">(<a href="#top">back to top</a>)</p>



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


<!-- CONTACT -->
## Contact

Julián Montoto-Louzao - [@julimontoto](https://twitter.com/julimontoto) - juli.mlouzao@gmail.com

GUANIN: [https://github.com/julimontoto/guanin](https://github.com/julimontoto/guanin)

<p align="right">(<a href="#top">back to top</a>)</p>




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png