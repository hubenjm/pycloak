from distutils.core import setup

DISTNAME = 'pycloak'
LICENSE = 'BSD'
MAINTAINER = "Mark Hubenthal"
EMAIL = "markhubenthal@gmail.com"
URL = "https://github.com/hubenjm/pycloak"
DESCRIPTION = (
"A framework and toolset for computing the double layer potential densities to "
"cancel the incident field on a control region. Restricted to the context of "
"the Helmholtz equation for now...")
PACKAGES = ['pycloak', 'plot']
setup(
name=DISTNAME,
version='1.0',
packages=PACKAGES,
license=LICENSE,
url=URL,
maintainer_email=EMAIL,
maintainer=MAINTAINER,
description=DESCRIPTION,
long_description=open('README.md').read()
)

