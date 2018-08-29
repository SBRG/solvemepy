from os import path

from setuptools.command.develop import develop
from numpy.distutils.core import Extension, setup

compile_args_quad = {"libraries": ["quadminos"]}
compile_args_double = {"libraries": ["minos"]}

if path.isfile("libquadminos.a"):
    compile_args_quad["library_dirs"] = [path.abspath(".")]
else:
    raise Exception('Missing libquadminos.a')

if path.isfile("libminos.a"):
    compile_args_double["library_dirs"] = [path.abspath(".")]
else:
    raise Exception('Missing libminos.a')

ext_modules = [
    Extension(name="qminospy.qwarmLP",
              sources=["qminospy/src/lp/qwarmLP.f90"], **compile_args_quad),
    Extension(name="qminospy.warmLP",
              sources=["qminospy/src/lp/warmLP.f90"], **compile_args_double),
    Extension(name="qminospy.qvaryME",
              sources=["qminospy/src/fva/qvaryME.f90"], **compile_args_quad),
    Extension(name="qminospy.qsolveME",
              sources=["qminospy/src/nlp/qsolveME.f90",
                       "qminospy/src/nlp/qmatrixA.f90"], **compile_args_quad),
    ]

setup(
    name="qminospy",
    ext_modules=ext_modules,
    cmdclass={"develop": develop}
    )
