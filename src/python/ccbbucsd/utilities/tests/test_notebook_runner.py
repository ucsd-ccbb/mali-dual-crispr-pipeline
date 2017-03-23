# standard libraries
import os
import random
import tempfile
import unittest

# test library
import ccbbucsd.utilities.notebook_runner as ns_test


class TestFunctions(unittest.TestCase):
    # No tests for read_in_notebook or get_params_from_nb_file as they just combine simple file opening and
    # calls to nbformat and nbparameterise; I'd essentially be writing tests for those third-party modules

    def get_html_subset(self, subset_name):
        orig_subset_1 = """<body>
  <div tabindex="-1" id="notebook" class="border-box-sizing">
    <div class="container" id="notebook-container">

<div class="cell border-box-sizing text_cell rendered">
<div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<p>A non-code cell</p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[&nbsp;]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">x</span> <span class="o">=</span> <span class="mi">5</span>
<span class="n">y</span> <span class="o">=</span> <span class="s2">&quot;blue&quot;</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[&nbsp;]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="nb">print</span><span class="p">(</span><span class="s2">&quot;x is </span><span class="si">{0}</span><span class="s2"> and y is </span><span class="si">{1}</span><span class="s2">&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="n">y</span><span class="p">))</span>
</pre></div>

</div>
</div>
</div>

</div>
    </div>
  </div>
</body>"""

        updated_subset_1 = """<body>
  <div tabindex="-1" id="notebook" class="border-box-sizing">
    <div class="container" id="notebook-container">

<div class="cell border-box-sizing text_cell rendered">
<div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<p>A non-code cell</p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[1]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">x</span> <span class="o">=</span> <span class="mi">6</span>
<span class="n">y</span> <span class="o">=</span> <span class="s1">&#39;blue&#39;</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[2]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="nb">print</span><span class="p">(</span><span class="s2">&quot;x is </span><span class="si">{0}</span><span class="s2"> and y is </span><span class="si">{1}</span><span class="s2">&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="n">y</span><span class="p">))</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area"><div class="prompt"></div>
<div class="output_subarea output_stream output_stdout output_text">
<pre>x is 6 and y is blue
</pre>
</div>
</div>

</div>
</div>

</div>
    </div>
  </div>
</body>"""

        updated_subset_2 = """<body>
  <div tabindex="-1" id="notebook" class="border-box-sizing">
    <div class="container" id="notebook-container">

<div class="cell border-box-sizing text_cell rendered">
<div class="prompt input_prompt">
</div>
<div class="inner_cell">
<div class="text_cell_render border-box-sizing rendered_html">
<p>A non-code cell</p>

</div>
</div>
</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[1]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="n">x</span> <span class="o">=</span> <span class="mi">5</span>
<span class="n">g_run_prefix</span> <span class="o">=</span> <span class="s1">&#39;faked_run_prefix&#39;</span>
</pre></div>

</div>
</div>
</div>

</div>
<div class="cell border-box-sizing code_cell rendered">
<div class="input">
<div class="prompt input_prompt">In&nbsp;[2]:</div>
<div class="inner_cell">
    <div class="input_area">
<div class=" highlight hl-ipython3"><pre><span></span><span class="nb">print</span><span class="p">(</span><span class="s2">&quot;x is </span><span class="si">{0}</span><span class="s2"> and g_run_prefix is </span><span class="si">{1}</span><span class="s2">&quot;</span><span class="o">.</span><span class="n">format</span><span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="n">g_run_prefix</span><span class="p">))</span>
</pre></div>

</div>
</div>
</div>

<div class="output_wrapper">
<div class="output">


<div class="output_area"><div class="prompt"></div>
<div class="output_subarea output_stream output_stdout output_text">
<pre>x is 5 and g_run_prefix is faked_run_prefix
</pre>
</div>
</div>

</div>
</div>

</div>
    </div>
  </div>
</body>
</html>"""

        if subset_name == "orig_1":
            result = orig_subset_1
        elif subset_name == "updated_1":
            result = updated_subset_1
        elif subset_name == "updated_2":
            result = updated_subset_2
        else:
            raise ValueError("Unrecognized subset_name '{0}'".format(subset_name))

        return result

    def get_temp_nb_str(self, get_str_1):
        nb_str_1 = """{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "A non-code cell"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "x = 5\\n",
    "y = \\"blue\\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": false
   },
   "outputs": [],
   "source": [
    "print(\\"x is {0} and y is {1}\\".format(x, y))"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.4.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 1
}
"""
        nb_str_2 = """{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "A non-code cell"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": true
   },
   "outputs": [],
   "source": [
    "x = 5\\n",
    "g_run_prefix = \\"blue\\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "collapsed": false
   },
   "outputs": [],
   "source": [
    "print(\\"x is {0} and g_run_prefix is {1}\\".format(x, g_run_prefix))"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.4.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 1
}
"""
        return nb_str_1 if get_str_1 else nb_str_2

    def write_temp_nb_file(self, get_str_1=True, parent_dir=None):
        input_nb_obj = tempfile.NamedTemporaryFile(dir=parent_dir, suffix=".ipynb")
        return self.write_nb_str_to_file_obj(input_nb_obj, get_str_1)

    def write_nb_str_to_file_obj(self, input_nb_obj, get_str_1=True):
        input_nb_str = self.get_temp_nb_str(get_str_1)
        input_nb_obj.write(bytes(input_nb_str, 'utf-8'))
        input_nb_obj.seek(0)
        return input_nb_obj

    def test_set_parameters(self):
        input_dict = {"x": 6,
                      "z": "red"}

        expected_output = 'x = 6\ny = \'blue\''

        input_nb_obj = self.write_temp_nb_file()
        # Note: not testing read_in_notebook here, just using it
        nb = ns_test.read_in_notebook(input_nb_obj.name)
        real_output = ns_test.set_parameters(nb, input_dict)

        # Ensure that: parameters in both notebook and input dictionary are replaced;
        # parameters not in input dictionary but in notebook are left as-is;
        # parameters in input dictionary but not in notebook are ignored
        self.assertEqual(expected_output, real_output.cells[1].source)

    # NOTE: for some reason that is unclear to me, this method *errors*--unable to start the kernel--when
    # run through the debugger, but *works* when run without the debugger.  Apparently kicking off a notebook
    # run with an outside process (nbconvert.ExecutePreprocessor) is ok but kicking off a kick off of a notebook run
    # with an outside process with another outside process (the debugger) is one layer of crazy too far.
    def test_execute_notebook(self):
        input_dict = {"x": 6,
                      "z": "red"}

        input_nb_obj = self.write_temp_nb_file()
        input_nb_dir, input_nb_filename = os.path.split(input_nb_obj.name)
        output_nb_filename = "{0}.ipynb".format(random.randint(1000000, 9999999))
        output_nb_fp = os.path.join(input_nb_dir, output_nb_filename)

        try:
            ns_test.execute_notebook(input_nb_filename, input_dict, output_nb_fp, run_path=input_nb_dir)

            # Note: not testing read_in_notebook here, just using it
            output_nb = ns_test.read_in_notebook(output_nb_fp)
        finally:
            if os.path.exists(output_nb_fp):
                os.remove(output_nb_fp)

        # test that the new variable has been correctly inserted and used to generate output in subsequent code cells
        self.assertEqual("x is 6 and y is blue\n", output_nb.cells[2].outputs[0]["text"])

    def test_export_notebook_to_html(self):
        input_nb_obj = self.write_temp_nb_file()
        input_nb_dir, input_nb_filename = os.path.split(input_nb_obj.name)

        output_fp = None
        try:
            output_fp = ns_test.export_notebook_to_html(input_nb_obj.name, input_nb_dir)

            with open(output_fp) as f:
                real_output = f.read()

        finally:
            if os.path.exists(output_fp):
                os.remove(output_fp)

        # Note: I don't check the entire HTML file contents, because it is huge and most of it is IPython-generated
        # CSS that (a) I don't care about and (b) might change.  I only check that the HTML I care about is a subset
        # of what is written to the HTML file.
        self.assertTrue(self.get_html_subset("orig_1") in real_output)


