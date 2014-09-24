Polyaffine Transformation Trees
===============================

For background information see [this](http://www.stanford.edu/~cseiler/).

To build binaries from the source it is required to install the latest ITK version.<br><br>

The code is based on the <a href="http://dx.doi.org/10.1016/j.neuroimage.2008.10.040">Log-Demons algorithm</a> (available on the <a href="http://hdl.handle.net/10380/3060">insight journal website</a> but not required here).

The main changes that we made are isolated to one method in one class,

<pre>itk::PolyaffineLogDomainDeformableRegistrationFilter::SmoothGivenField(VelocityFieldType * field, <br>const double StandardDeviations[ImageDimension]) </pre>

This function is called to regularize the estimated SVF after every iteration of the log-demons algorithm. We overwrite this method to estimate the polyaffine transformation tree. The first argument is the unregularized SVF. The second argument is the amount of smoothing applied in the standard Gaussian smoothing step. In our implementation the second arguments is not used.<br><br>

<b>Examples:</b><br>

The following example command registers a template image to a subject with five levels (-s 5) and 2 voxel dilation of the union of the mask images (-r 2),

<pre>> PolyaffineLogDomainDemonsRegistration -f TemplateImage.mhd -F TemplateMaskImage.mhd
-m SubjectImage.mhd -M SubjectMaskImage.mhd -s 5 -r 2</pre>

Registration with prior mean (--prior-Mean PriorMean.txt) and concentration matrix (--prior-Concentration PriorConcentration.txt) on the tranformation parameter, and 4 mm variance of the velocity vectors (&#8209;w 4),

<pre>> PolyaffineLogDomainDemonsRegistration -f TemplateImage.mhd -F TemplateMaskImage.mhd
-m SubjectImage.mhd -M SubjectMaskImage.mhd -s 5 -r 2
--prior-Mean PriorMean.txt --prior-Concentration PriorConcentration.txt -w 4</pre>

Registration with manual regions (-R LabelImageManualRegions.mhd, the labels are expected to be numbered with integers ranging from one to the total number of regions),

<pre>> PolyaffineLogDomainDemonsRegistration -f TemplateImage.mhd -F TemplateMaskImage.mhd
-m SubjectImage.mhd -M SubjectMaskImage.mhd -s 2 -r 2 -R LabelImageManualRegions.mhd</pre>

In case of questions, please don't hesitate to contact me.
