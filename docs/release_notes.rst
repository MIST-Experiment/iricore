Release notes
-------------

v1.8.2
======
* Fixed binary compilation bug during installation on MacOS

v1.8.1
======
* Implemented autoupdate of indices
* Fixed reading IRI data to memory in master Python process - now `iricore`
  should be faster when called repeatedly in a loop.
* Small bugfixes and refactoring.

v1.8.0
======
* Updated IRI-2020 to the 03/05/2024 version.
* Added support for OARR IRI output. It is now included in the returned :class:`iricore.IRIOutput`.
* Added support for the manual user input for IRI parameters - see `**kwargs` in :func:`iricore.iri`.

v1.7.0
======
* Added :func:`iricore.refstec` that includes simple model of ray refraction in the slant TEC calculation.
* :func:`iricore.stec` and :func:`iricore.refstec` now allow a custom array of heights as an input. They now also can
  return a ray trajectory and additional info if `return_hist` parameter is True.


