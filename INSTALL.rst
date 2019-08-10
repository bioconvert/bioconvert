Note for developers.

biocode is not installable from conda so this is causing trouble with travis.
since it is used in one converter only, we removed it from the installation.
Yet, it should be installed by developer. Hence in singularirty and Docker file,
we install biocode manually
