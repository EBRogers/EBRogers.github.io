---
title: 'Postmortem: Adding contrast testing to Cardinal'
date: 2025-10-25
permalink: /posts/contrasts
tags: Cardinal MSI
---

{% include toc %}

[Bioconductor 3.22](https://bioconductor.org/developers/release-schedule/) is being released on October 21, 2025. Among bug fixes, this release will include some useful additions to [Cardinal for MSI](https://cardinalmsi.org), specifically the addition of contrast testing for linear models. I [added](https://github.com/kuwisdelu/Cardinal/pull/69) this and wanted to detail some of my thoughts and decisions below.

<details markdown="1">
<summary><strong>What is this for? What is a contrast?</strong></summary>
Some real quick background: I do statistical analysis of mass spectrometry imaging (MSI) data. Basically the stuff after the data are collected, preprocessed, and ROIs selected - like *modeling* and *hypothesis testing*. In order to make *population-level statements* from a sample of MSI images (e.g., "on average, exposure to this drug increases the abundance of this analyte") you have to specify the relationships using a statistical model and use them to estimate the signal and variance. MSI data have many sources of variance so we need multiple samples (biological replicates). Also to make statements about the population effect from our sample we need models that specify what we think is happening (linear models), a way to test the signal we observe (the effect) against the noise of the effect (the variance), and an empirical way to compare that signal-to-noise ratio (a t statistic) to what we would expect to see based on random chance alone (a p value). You have probably done this - it's a t test<sup>1</sup>.

A contrast specifies a specific comparison from a statistical model using the model parameters. Let's go back to the t test as an example. Take a two sample t test comparing the average of treatment vs control.

$$
\frac{\bar{Y}_{Treatment}-\bar{Y}_{Control}}{\sqrt{\frac{2\sigma^2}{n}}} \sim t_{df}
$$

This test actually assumes a linear model beneath it - or a one-way ANOVA (it's all the same math).

$$
Y_i = \beta_0 + X_i\beta_1 + \epsilon_i \text{ where } \epsilon_i \overset{iid}{\sim} N(0,\sigma^2)
$$

Rewritten in a more friendly way:

$$
\displaylines{
Y_i = \text{Mean Y for Control} + X_i(\text{Mean difference of Y between control and treatment}) + \text{Error}_i \\\
\text{where } X_i = 0 \text{ for control, or } 1 \text{ for treatment}
}
$$

A contrast is a set of coefficients that usually sum to 0, and correspond to model parameters. In abstract, it's how we go from the fitted model to a specific estimate. For the t test above, comparing the difference between treatment and control, the contrast would be \\((0,1)\\) and the estimate being \\((0\times\beta_0 + 1\times\beta_1)\\) where \\(\beta_0\\) is the average response of the control and \\(\beta_1\\) is the average difference between the control and treatment, both from the fitted model. Follow me?

Now for a simple example, the two sample t test and the contrast with a linear model are the same. The real magic happens when you have multiple independent variables, levels of each variable, sources of variation (like random factors in a mixed-effects model), and interactions. The contrast also specifies how to compute the *variance* of the estimate. I've heard this referred to as "propagating error" for people familiar with GraphPad Prism, and, in essence, it lets you piggyback off of the explanatory power of the other variables, random effects, and interactions to get a more precise estimate of variance.


<sup>1</sup>In the purest form, a t test is a signal-to-noise ratio. A special case is a two sample t test or a paired t test where the signal is the difference of two means and the ratio is the pooled variance.

</details>

## What did Cardinal for MSI have before?

`Cardinal::meansTest()` did and does do differential analysis, depending on exactly what difference you are trying to analyze. `meansTest()` computes mean intensity for each feature and sample then fits linear models, with `stats::lm()` or `nlme::lme()` when a random effect is specified. `meansTest()` then performs a *likelihood ratio test* (LRT) against a reduced model, defaulting to `~ 1`. The interpretation of this test changes a lot depending on the context, but in general an LRT compares a full and reduced model to test if the full model explains significantly more variation than the reduced. In this use case, an LRT is equivalent to an F test - and if you are testing a model with only a single fixed effect with two levels, it's equivalent to a two sample t test. But beyond that situation (testing a single fixed effect) LRTs quickly become less useful if your goal is to test differential abundance. For one thing, you don't get estimates of the effect - just if it's significant. Also you can't test individual factor levels.

So instead let's t test (signal to noise ratio remember?) contrasts of interest. Now how to do that... you could, of course, extract the fitted models from the `MeansTest` result object and run the tests yourself. This is fine for models without random effects that use `lm` but `nlme::lme` models don't check for singularity (in this case, most often when subject/biological variation is too small to estimate) and offer a [primitive method of degrees of freedom computation](https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#why-doesnt-lme4-display-denominator-degrees-of-freedomp-values-what-other-options-do-i-have). Is that a big deal? Maybe - probably not if you're doing something simple - but might as well compute means within Cardinal and use `lme4::lmer()` on the `FeatureData` frame pulled from every sample, seeing as how you already have to run the contrasts yourself anyway. `lme4::lmer()` lets you know when the models are singular so you can check for (then ignore :sunglasses:) zero subject-level variance and `lmerTest` offers a method for Satterthwaite approximation of the degrees of freedom. Contrast testing can be done with `emmeans`. It's complex, but if nothing else, `emmeans`'s wide adoption yields many examples of its usage with different research questions and experimental designs.

But in my experience, most Cardinal users are experimentalists, not computer scientists or statisticians, and won't or can't do the extra steps above because they are too complicated. It makes sense - we all have our areas of expertise and no one can be specialized in everything. So since I have this expertise, and the need is clear (at least to me), I decided to add contrast testing to Cardinal.

## Cardinal for MSI does contrasts now
Running contrasts in Cardinal is relatively straightforward. All MSI samples need to be in the same `MSImagingExperiment` object, and pixels from each sample identified using a factor in the `PositionDataFrame`, usually the run. Independent variables must be specified in the same way. Then it's pretty simple - first, run `meansTest` with your preferred specification setting `use_lmer = TRUE`. Then feed the results into `contrastTest`, specifying the contrasts using `specs` and `method` args. These pretty much map directly to `specs` arg of `emmeans::emmeans()` and the `method` arg of `emmeans::contrast`. Exactly how to specify these depends on your hypotheses of interest, but as I mentioned above emmeans has extensive [documentation](https://rvlenth.github.io/emmeans/), [vignettes](https://cran.r-project.org/web/packages/emmeans/index.html), and [examples](https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/). Or just ask ChatGPT, with appropriate skepticism. `Cardinal::topFeatures` will do FDR correction, the very important final step. I'm finishing up a paper about this, available very soon, that will include code and workflows for you to use with your data.

## Thoughts and caveats on the implementation

Although I am a relatively competent scientific programmer (I write code to do science pretty well) that is not the same as being a good software engineer... We all hate on vibe coding (god I hate that word) and all its cringy etymological children (search LinkedIn for "vibe coding cleanup specialist" to see truly the most uninspired people), but in their narrow and well confined lane I found Windsurf, Cursor, and GitHub Copilot to be very useful. It helped me plan, answered questions about Cardinal and its existing repository, and frankly "thought" of stuff I would have missed. My expertise actually came with testing - I know how these models work and what the results should be. I have used Cardinal a lot and I know how I want it to behave. I can create realistic test data for each scenario and I had already implemented the non-Cardinal version of this process from above. 

Now onto some interesting technical details in no particular order. Upon conference with Kylie (my colleague and the creator of Cardinal), we decided to keep `nlme::lme()` and simply add `lme4::lmer()`, specified with the `use_lmer` arg. Why? First we wanted the LRT to behave exactly the same as it always has. `lme4` and `nlme` are relatively interchangeable with mostly esoteric differences, but their estimators are not exactly the same so they don't produce the *exact* same answers. Better to keep the old way around in case someone out there relies on LRTs. Second, `nlme` supports some interesting additions, such as [Matérn covariance structures](https://search.r-project.org/CRAN/refmans/spaMM/html/MaternCorr.html) that could let us specify pixel-level models that overcome the biasing effect of [pixel autocorrelation](https://doi.org/10.1021/acs.analchem.6b00672) inherent to MSI data. In the same vein, I decided not to implement LRT testing with `lme4::lmer()` models. `lme4` offers both maximum likelihood estimation (MLE) and restricted maximum likelihood estimation (REML - I know but it's convention). MLE generally underestimates variance, especially at small sample sizes whereas REML is generally considered unbiased ([see this](https://doi.org/10.1080/00273171.2017.1344538)). Unfortunately, LRTs with REML-estimated parameters must have the same fixed effects, rendering them useless for our purposes. As MSI experiments often have few (if any) replicates, in my opinion, it's more important to use REML than it is to offer LRTs.

## Conclusions
Going forward, I think my contributions to Cardinal will focus on filling usability gaps. Contrast testing seemed to be the largest hole in the end to end workflow of differential analysis - but there are others, mostly small and avoidable, but a few (deisotoping) that still need to be remedied with external packages. Excited for the future.