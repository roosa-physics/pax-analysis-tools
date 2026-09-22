Adding a detector response model
--

After getting an input spectrum we need to apply a detector response model.

I approximate this with a gaussian but we have a problem. If I calculate a full convolution in a fitting routine, the sun will die before my fit completes.

We solve this by making a look-up table.