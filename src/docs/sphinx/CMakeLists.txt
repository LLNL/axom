################################
# documentation targets
################################

set(SPHINX_HTML_THEME "
try:
    import sphinx_rtd_theme
except:
    html_theme = 'classic'
    html_theme_options = {
        'codebgcolor': 'lightgrey',
        'stickysidebar': 'true'
    }
    html_theme_path = []
else:
    html_theme = 'sphinx_rtd_theme'
    html_theme_options = {}
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
")

add_subdirectory(web_main)
add_subdirectory(coding_guide)
add_subdirectory(dev_guide)
add_subdirectory(quickstart_guide)



