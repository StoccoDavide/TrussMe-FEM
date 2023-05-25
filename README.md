# TrussMe FEM

A Maple package for symbolic FEM structural analysis.

## Installation

To install the module you must have first installed Maple. Then open the `PackAndGo.mw` file and use the `!!!` button to *execute the entire worksheet*.

Then test the module in a Maple worksheet or document by typing:

```
> TrussMe_FEM:-Info(TrussMe_FEM); # For Maple versions up to 2020
> TrussMe_FEM:-Info();    # For Maple versions starting from 2021
```

Alternatively, you can use one of the tests file provided in the `tests` folder. If the module is loaded without errors, it is done!

## Module description

If you want a full description of the `TrussMe_FEM` module:

```
> Describe(TrussMe_FEM);
```

This command will generate a brief description of the module and all the procedures and other objects present in `TrussMe_FEM.mpl` file.

## Usage

In case you have no time to read the description and realize how it should or should not work, refer to the files in the `tests` folder.

### 🚧 Attention! 🚧

Maple object-oriented programming features have slightly changed in 2021, which online documentation states:

> As of Maple 2021, if the method has a formal parameter named `_self`, references to its object's local or exported variables may be written without prefixing them. That is, `_self:-variable` may be written as just `variable`. Maple will add the `self:-` prefix internally when the method is simplified.

> As of Maple 2021, a message-passing form of method call can be used, which will automatically pass the object as an argument if the method has a formal parameter named `_self`.

> Another way to invoke a method, similar to that used in other object-oriented languages, was introduced in Maple 2021. In this form, the object name is qualified by the method name, `object_name:-method_name(argument)` just as it can be in the function mechanism described above. However, the *object can be omitted from the argument sequence*.

For further information please refer to this [link](https://fr.maplesoft.com/support/help/Maple/view.aspx?path=object/methods).

## Authors

- *Davide Stocco*,
  Department of Industrial Engineering,
  University of Trento \
  email: davide.stocco@unitn.it

- *Matteo Larcher*,
  Department of Industrial Engineering,
  University of Trento \
  email: matteo.larcher@unitn.it
