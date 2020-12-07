import pathlib
from distutils.core import setup

HERE = pathlib.Path(__file__)

README = (HERE / "README.md").read_text()

setup(
    name="abxrxpro",
    version="2.0.1",
    description="Visualise your isolate's antibiotic resistance profile",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/CaileanCarter/AbxRxPro",
    author="Cailean Carter",
    author_email="cailean.carter@quadram.ac.uk",
    license="MIT",
    include_package_data=True,
    install_requires=["plotly", "pandas"] 
)