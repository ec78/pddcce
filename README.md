# GAUSS Panel Data Dynamic Common-Correlated Effects Model
Econometric package for performing panel data dynamic common-correled effects model (DCCE).

## Getting Started
### Prerequisites
The DCCE library requires:

1.  A working copy of **GAUSS 23+**.


### Installing
The GAUSS DCCE library can be installed using the [GAUSS applications installer](https://www.aptech.com/support/installation/using-the-applications-installer-wizard/).

**Please do not download the source code and install manually. You will not be able to properly install the library.**

Before using the functions created by `dccelib` you will need to load the newly created `dccelib` library. This can be done in a number of ways:
  *  Navigate to the library tool view window and click the small wrench located next to the `dccelib` library. Select `Load Library`.  
  *  Enter `library dccelib` in the program input/output window.
  *  Put the line `library dccelib;` at the beginning of your program files.

>  Note: I have provided the individual files found in `sslib_1.0.zip` for examination and review. However, installation should always be done using the GAUSS Package Manager.

### Examples
After installing the library, examples for all available procedures can be found in your **GAUSS** home directory in the directory **pkgs > dccelib >examples**.

The following example files are currently included in the GAUSS dccelib:

|File      |     Description              |
|----------|------------------------------|
| mg_penn.e | Example mean group estimation using a panel dataset from the Penn World Tables. The data ranges from 1960 to 2007.|
|cce_penn.e | Example common-correlated effects estimation using a panel dataset from the Penn World Tables. The data ranges from 1960 to 2007.|
|dcce_penn.e | Example dynamic common-correlated effects estimation using a panel dataset from the Penn World Tables. The data ranges from 1960 to 2007.|

## License
The author makes no performance guarantees. The `dccelib` is available for public non-commercial use only.

## Author
For any bugs, please send e-mail to [Eric Clower](mailto:eric@aptech.com).

[Aptech Systems, Inc](https://www.aptech.com/)  
[![alt text][1.1]][1]
[![alt text][2.1]][2]
[![alt text][3.1]][3]

<!-- links to social media icons -->
[1.1]: https://www.aptech.com/wp-content/uploads/2019/02/fb.png (Visit Aptech Facebook)
[2.1]: https://www.aptech.com/wp-content/uploads/2019/02/gh.png (Aptech Github)
[3.1]: https://www.aptech.com/wp-content/uploads/2019/02/li.png (Find us on LinkedIn)

<!-- links to your social media accounts -->
[1]: https://www.facebook.com/GAUSSAptech/
[2]: https://github.com/aptech
[3]: https://linkedin.com/in/ericaclower
<!-- Please don't remove this: Grab your social icons from https://github.com/carlsednaoui/gitsocial -->
