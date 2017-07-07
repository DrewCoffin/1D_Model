IO TORUS CHEMISTRY MODEL
(fixed-width readme file, see aareadme.pdf for a detailed readme)

CONTENTS
             model/: Dir containing the onebox model and JHPL IDL library
            widget/: Dir containing the widget, see readme file within widget dir.
        compile.pro: IDL script file that compiles the model, supporting functions, as well as defines system variables
       runmodel.pro: IDL script file that runs one instance of the model
dataset_creator.pro: IDL procedure that runs model through a wide range of parameter space


USAGE
To compile the model:
IDL> @compile

To run one instance of the model:
IDL> @compile
IDL> @runmodel

To create data set, modify parameters within procedure, then:
IDL> @compile
IDL> dataset_creator

