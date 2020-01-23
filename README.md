# DNP NET

[![N|Solid](http://www.bridge12.com/wp-content/uploads/2016/10/b12logo.png)](http://www.bridge12.com/)

DNP net is an open-source, python library for importing and processing ODNP data.

  - Type some Markdown on the left
  - See HTML in the right
  - Magic

# Features

  - Import ODNP Data from Bruker, Varian, and Kea formats
  - Construct n-dimensional arrays easily
  - Process Data including window, zero-filling, Fourier transforms, etc.

```python
data.window('t')
data.zero_fill('t',n_pts = 1024)
data.ft('t')
data.plot('t')
show()
```

You can also:
  - Export data in ASCII Format
  - Create quality figures
