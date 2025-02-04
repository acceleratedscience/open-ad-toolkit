import os
import base64
from IPython.core.display import HTML


class NotebookStyles:
    """Generate the required CSS to style OpenAD Notebooks"""

    dir_path = None  # Path of the parent directory
    css_main = "AAA"  # Holds the main CSS code
    css_banner_main = ""  # CSS code for main bannmer
    css_banner_chapter = ""  # CSS code for secondary banner
    css_output = ""  # Combined CSS for output

    def __init__(self):
        self.dir_path = os.path.dirname(os.path.abspath(__file__))
        self.load_main_css()
        self.create_banners_css()
        css = "\n".join([self.css_main, self.css_banner_main, self.css_banner_chapter])
        self.css_output = HTML(f"<style>{css}</style>")

    def css(self):
        """Return the final CSS code"""
        return self.css_output

    def load_main_css(self):
        """Load the CSS code from the main CSS file"""
        path_css = os.path.join(self.dir_path, "custom.css")
        self.css_main = open(path_css, "r", encoding="utf-8").read()

    def create_banners_css(self):
        """Generate CSS code for the banners"""
        path_banner_main = os.path.join(self.dir_path, "media", "banner-bg.jpg")
        path_banner_chapter = os.path.join(self.dir_path, "media", "science_banner.jpg")
        data_uri_banner_main = self._image_to_data_uri(path_banner_main)
        data_uri_banner_chapter = self._image_to_data_uri(path_banner_chapter)
        self.css_banner_main = self._create_banner_css("banner:not(.chapter)", data_uri_banner_main)
        self.css_banner_chapter = self._create_banner_css("banner.chapter", data_uri_banner_chapter)

    def _image_to_data_uri(self, filepath):
        """Convert an image to a data URI (url in css doesn't work)"""
        with open(filepath, "rb") as image_file:
            base64_encoded_data = base64.b64encode(image_file.read()).decode("utf-8")
            return f"data:image/jpeg;base64,{base64_encoded_data}"

    def _create_banner_css(self, class_name, data_uri):
        return (
            f".{class_name}::after {{ background: url({data_uri}) center center no-repeat; background-size: cover; }}"
        )
