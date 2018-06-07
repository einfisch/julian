#Example_a3: Prolate spheroidal coordinates, fibres: viewing the heart
#
# The heart model uses prolate spheroidal coordinates to describe its complex,
# rounded shape with minimal degrees of freedom. The mesh still requires many
# elements, and the heart wall is two elements thick in places. In addition,
# fibre fields describe the orientation of muscle fibres and sheets over the
# heart. This example creates graphics to help visualise the structure of the
# heart, demonstrates prolate spheroidal coordinates, and shows how to view
# a fibre field.
#
# Be sure to run the commands in this example in the order they are presented,
# and not all at once. Comments and instructions for things you can do are given
# between these commands.
#
# The exnode and exelem files used in this example were generated from 
# example 141
#
#----------
#
# Create a few materials in addition to the default for use later.
gfx create material heart ambient 0.3 0 0.3 diffuse 1 0 0 specular 0.5 0.5 0.5 shininess 0.5;
gfx cre mat trans_purple ambient 0.4 0 0.9 diffuse 0.4 0 0.9 alpha 0.3
gfx cre mat bluey ambient 0 0.2 0.4 diffuse 0 0.5 1 specular 0.5 0.5 0.5 shininess 0.8
gfx cre mat gold ambient 1 0.7 0 diffuse 1 0.7 0 specular 1 1 0.8 shininess 0.8
gfx cre mat axes ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5
#
# Read in the heart and draw axes.
gfx read nodes example heart.exnode
gfx read elements example heart.exelem
gfx create axes material axes length 1.5
gfx draw axes
#
# Open the graphics window and reorient the heart up the right way (since -x is
# "up" in the heart model). Also turn on perspective.
gfx cre win 1
gfx mod win 1 image rotate 0 1 0 -90
gfx mod win 1 view perspective
#
# Draw the endo-cardial surfaces in the heart colour. More exactly, all exterior
# surfaces of the model for which coordinate xi3 = 0.
gfx mod g_e heart surfaces exterior face xi3_0 mat heart
#
# Now we'll spend some time looking at iso-surfaces of the prolate spheroidal
# coordinates. Firstly, define a field for each of the components of the coordinates: lambda, mu, theta
gfx define field coordinates.lambda component coordinates.lambda
gfx define field coordinates.mu component coordinates.mu
gfx define field coordinates.theta component coordinates.theta
#
# Next open the scene editor. 
gfx edit scene
#
# Set 'Auto' to off, so that 'Apply' will be needed to see the iso-surfaces.
# The group 'heart' on scene default should be selected, and it will have the
# default lines and the endo-cardial surfaces added earlier. Select
# 'iso_surfaces', 'Add' and make sure the iso-scalar 'coordinates.lambda' is to
# equal to 0.5 before clicking 'Apply'. Coordinate lambda increases in curves
# away from the central axis of the heart, (-)x. Change the iso value to 0.75
# and apply the changes. Next change the iso-component to mu, which increase up
# the heart. Set mu to 1.5708 (=PI/2) to see a nice slice through the heart.
# Finally, try the theta component which varies from 0 to 2*PI around the axis
# of the heart.
#
# Try having several different iso-surfaces visible simultaneously, and with
# different materials. Change the element discretization (in the general
# settings) to 8*8*4 and apply this. It will take a few seconds to compute the
# graphics at this higher quality rendition, but the surfaces and iso-surfaces
# should look somewhat nicer. Note that due to the shape of the mesh and the
# elements used, it is fine to use less detail to draw graphics in the xi3
# (radial/lambda) direction.
#
# Now add semi-transparent purple surfaces to the exterior of the heart to see
# its outer shape and interior simultaneously.
gfx mod g_e heart surfaces exterior face xi3_1 mat trans_purple
# (If you haven't enabled the "Auto Apply/Revert" feature of the scene editor,
# you'll have to click on "Revert" to reload the rendition from the scene into
# the editor.)
#
# For the next part of this example. select all the iso-surfaces one-by-one and
# click on the 'Del' key to delete them. Also make the outside surfaces
# invisible, then click 'Apply'.
#
# Add element_points to the rendition with the discretization set at 4*4*1 and
# the material gold. Choose the 'cylinder' glyph, make the base glyph size
# 2*0.5 (becomes 2*0.5*0.5 when you press enter), choose field 'fibres' as the
# orientation/scale field and set the scale_factors (multiplying the magnitude
# of the orientation/scale field vectors) to 0, since the fibre axes are unit
# vectors specifying direction only. The heart is now covered in small cylinders
# indicating the direction of the fibres. Following the discretization values,
# there are 4*4*1 cylinders across each element in its xi1*xi2*xi3
# directions, respectively. Note that the 'Glyph_size' paramaters display the
# size as length*width*height - the width is in the direction of the fibre
# sheet.
#
# Cylinders are just one way of viewing the fibres; they don't show the plane of
# the sheet, nor do they indicate which direction they are pointing in along
# their axis - which may be useful for checking your fibre field. You can choose
# several other "Glyphs" besides cylinders to display over the fibre field.
# Apply the 'arrow_solid' glyph and take a close look at it. (The centre of the
# cylinder is half-way along its length, whereas the arrow is centred at its
# tail. The 'Centre' can be changed on the editor. All glyphs are designed to be
# of unit dimension - 1 unit long, 1 unit diameter at the widest point etc. -
# so sensible centre values should be of a similar magnitude.)
#
# Go back to using cylinders, add another set of fibres using 'sheets' and
# material 'bluey'. Enter '1.5' in the 'Glyph size' field. You can now see both
# the the fibre direction and sheet plane. Experiment with different glyphs,
# sizes, discretizations and materials to see how you can interpret the fibre
# field. In general, keep the discretization the same for all fibre graphics.
# 
#----------
#
