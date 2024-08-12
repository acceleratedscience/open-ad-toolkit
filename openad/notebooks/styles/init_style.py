import base64
from IPython.core.display import HTML


def init():
    """
    Initialize some styling that can be reused across Notebooks.
    """

    # Load general CSS
    def css_styling():
        styles = open("./styles/custom.css", "r").read()
        return styles

    style1 = css_styling()

    ###
    ###

    # Load banner image as data (url doesn't work).
    def image_to_data_uri(filepath):
        """
        Convert an image to a data URI.
        """
        with open(filepath, "rb") as image_file:
            base64_encoded_data = base64.b64encode(image_file.read()).decode("utf-8")
            return f"data:image/jpeg;base64,{base64_encoded_data}"

    # Use this function to get the data URI for your image
    data_uri = image_to_data_uri("./media/banner-bg.jpg")

    # Then, in your CSS, you can use this data URI directly
    def banner_styling():
        styles2 = f"""
        <style>
        .banner::after {{
            background: url({data_uri}) center center no-repeat;
            background-size: cover;
        }}
        </style>
        """
        return styles2

    style2 = banner_styling()

    return HTML(style1 + style2)
