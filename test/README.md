# Polypolish tests

Polypolish comes with a few automated tests to help with development and spotting bugs.

To run the tests, execute this command from Polypolish's root directory:
```
pytest
```

Or if you have [Coverage.py](https://coverage.readthedocs.io/en/coverage-4.5.1a/) installed, you can run the tests through it to get some code coverage stats:
```
coverage run -m pytest && coverage report -m
```
