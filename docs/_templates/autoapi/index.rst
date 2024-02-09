API Reference
=============

This page contains auto-generated API reference documentation, which describes
pygscatalog libraries.

The information is mostly useful for developers: people that want to write Python
code to work with polygenic scores.


.. toctree::
   :titlesonly:

   {% for page in pages | sort %}
   {#
      Add the top most levels in "pgscatalog.X" to the index file
      This is needed because we don't have __init__.py file in pgscatalog package
      as we use nested implicit namespace packages.
      https://github.com/readthedocs/sphinx-autoapi/issues/298
   #}
   {% if (page.top_level_object or page.name.split('.') | length == 2) and page.display %}
   {{ page.short_name }} <{{ page.include_path }}>
   {% endif %}
   {% endfor %}
