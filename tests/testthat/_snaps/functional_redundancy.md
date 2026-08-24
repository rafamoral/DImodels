# functional redundancy works

    Code
      print(AV_res)
    Output
      -------------------------------------------------------------------------------- 
      
      Testing functional redundancy for species p3 and p4.
      
      -------------------------------------------------------------------------------- 
      
      Equivalence tests conducted at the alpha = 0.05 level (i.e, using 90% CIs) with a delta margin of (-1.27, 1.27) quantified using 0.75 times residual SD of model.
      
      Null hypothesis (H0): The species are not functionally redundant.
      Alternative hypothesis (H1): The species are functionally redundant.
      
      --------------------------------------------------------------------------------
      
      Results: 
      
      Functional redundancy was not established for p3 and p4 (p = 1.00) due to the following reasons.
      
        - Monoculture performances were not found to be equivalent (p = 0.99).
        - Within-group interactions were not found to be equivalent to 0 (p = 1.00).
        - Between-group interactions were not found to be equivalent for the following (p = 0.87):
          * p1 (p = 0.87)
          * p2 (p = 0.87)
      
      --------------------------------------------------------------------------------
      Note: All tests follow the intersection-union principle and the respective overall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      summary(AV_res)
    Output
      -------------------------------------------------------------------------------- 
      
      Summary of functional redundancy test for p3 and p4.
      
      -------------------------------------------------------------------------------- 
      
      Hypotheses:
      
        H0: The species are not functionally redundant.
        HA: The species are functionally redundant.
      
      Test specification:
      
        Alpha: 0.05
        Confidence level: 90%
        Equivalence margin: (-1.27, 1.27) [quantified using 0.75 times residual SD of model]
      
      Overall result:
      
        Functional redundancy not established (p = 1.00).
      
      Individual Criteria:
      
        Monoculture performances (p = 0.99): [0/1] contrasts equivalent to 0.
      
          Contrast  Estimate  p.value
          --------  --------  -------
          p3 - p4   3.485     0.99   
      
        Within-group interactions (p = 1.00): [0/1] contrasts equivalent to 0.
      
          Contrast        Estimate  p.value
          --------------  --------  -------
          equi - mono_av  3.568     1.00   
      
        Between-group interactions (p = 0.87): [0/2] contrasts equivalent to 0.
      
          Contrast       Estimate  p.value
          -------------  --------  -------
          p1_p3 - p1_p4  1.742     0.87   
          p2_p3 - p2_p4  1.742     0.87   
      
      Note:  All tests follow the intersection-union principle and the respectiveoverall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      print(FULL_res)
    Output
      -------------------------------------------------------------------------------- 
      
      Testing functional redundancy for species p1 and p2.
      
      -------------------------------------------------------------------------------- 
      
      Equivalence tests conducted at the alpha = 0.05 level (i.e, using 90% CIs) with a delta margin of (-5, 5) as defined by user.
      
      Null hypothesis (H0): The species are not functionally redundant.
      Alternative hypothesis (H1): The species are functionally redundant.
      
      --------------------------------------------------------------------------------
      
      Results: 
      
      The species p1 and p2 are functionally redundant (p = 0.0077).
      
      All equivalence criteria were satisfied.
      
        - Monocultures were found to have equivalent performance (p < 0.0001). 
        - Within-group interactions were found to be equivalent to 0 (p = 0.0077). 
        - Between-group interactions with all other species were found to be equivalent (p = 0.0004).
      
      --------------------------------------------------------------------------------
      Note: All tests follow the intersection-union principle and the respective overall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      summary(FULL_res)
    Output
      -------------------------------------------------------------------------------- 
      
      Summary of functional redundancy test for p1 and p2.
      
      -------------------------------------------------------------------------------- 
      
      Hypotheses:
      
        H0: The species are not functionally redundant.
        HA: The species are functionally redundant.
      
      Test specification:
      
        Alpha: 0.05
        Confidence level: 90%
        Equivalence margin: (-5, 5) [user specified delta]
      
      Overall result:
      
        Functional redundancy established (p = 0.0077).
      
      Individual Criteria:
      
        Monoculture performances (p < 0.0001): [1/1] contrast equivalent to 0.
      
          Contrast  Estimate  p.value
          --------  --------  -------
          p1 - p2   -0.249    0.0001 
      
        Within-group interactions (p = 0.0077): [1/1] contrast equivalent to 0.
      
          Contrast        Estimate  p.value
          --------------  --------  -------
          equi - mono_av  2.510     0.0077 
      
        Between-group interactions (p = 0.0004): [2/2] contrasts equivalent to 0.
      
          Contrast       Estimate  p.value
          -------------  --------  -------
          p3_p1 - p3_p2  0.316     0.0004 
          p4_p1 - p4_p2  0.183     0.0003 
      
      Note:  All tests follow the intersection-union principle and the respectiveoverall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      print(base_res)
    Output
      -------------------------------------------------------------------------------- 
      
      Testing functional redundancy for species p1 and p2.
      
      -------------------------------------------------------------------------------- 
      
      Equivalence tests conducted at the alpha = 0.05 level (i.e, using 90% CIs) with a delta margin of (-1, 1) as defined by user.
      
      Null hypothesis (H0): The species are not functionally redundant.
      Alternative hypothesis (H1): The species are functionally redundant.
      
      --------------------------------------------------------------------------------
      
      Results: 
      
      Functional redundancy was not established for p1 and p2 (p = 0.89) due to the following reasons.
      
        - Monoculture performances were not found to be equivalent (p = 0.40).
        - Within-group interactions were not found to be equivalent to 0 (p = 0.85).
        - Between-group interactions were not found to be equivalent for the following (p = 0.89):
          * p3 (p = 0.89)
          * p4 (p = 0.12)
      
      --------------------------------------------------------------------------------
      Note: All tests follow the intersection-union principle and the respective overall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      summary(base_res)
    Output
      -------------------------------------------------------------------------------- 
      
      Summary of functional redundancy test for p1 and p2.
      
      -------------------------------------------------------------------------------- 
      
      Hypotheses:
      
        H0: The species are not functionally redundant.
        HA: The species are functionally redundant.
      
      Test specification:
      
        Alpha: 0.05
        Confidence level: 90%
        Equivalence margin: (-1, 1) [user specified delta]
      
      Overall result:
      
        Functional redundancy not established (p = 0.89).
      
      Individual Criteria:
      
        Monoculture performances (p = 0.40): [0/1] contrasts equivalent to 0.
      
          Contrast  Estimate  p.value
          --------  --------  -------
          p1 - p2   0.743     0.40   
      
        Within-group interactions (p = 0.85): [0/1] contrasts equivalent to 0.
      
          Contrast        Estimate  p.value
          --------------  --------  -------
          equi - mono_av  -3.614    0.85   
      
        Between-group interactions (p = 0.89): [0/2] contrasts equivalent to 0.
      
          Contrast       Estimate  p.value
          -------------  --------  -------
          p3_p1 - p3_p2  -2.841    0.89   
          p4_p1 - p4_p2  0.372     0.12   
      
      Note:  All tests follow the intersection-union principle and the respectiveoverall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      print(nine_res_v1)
    Output
      -------------------------------------------------------------------------------- 
      
      Testing functional redundancy for species p1, p2, p3, p4, and p5.
      
      -------------------------------------------------------------------------------- 
      
      Equivalence tests conducted at the alpha = 0.05 level (i.e, using 90% CIs) with a delta margin of (-2, 2) as defined by user.
      
      Null hypothesis (H0): The species are not functionally redundant.
      Alternative hypothesis (H1): The species are functionally redundant.
      
      --------------------------------------------------------------------------------
      
      Results: 
      
      Functional redundancy was not established for p1, p2, p3, p4, and p5 (p = 1.00) due to the following reasons.
      
        - Monoculture performances were not found to be equivalent (p = 1.00).
        - Within-group interactions were not found to be equivalent to 0 (p = 0.31).
        - Between-group interactions were not found to be equivalent for the following (p = 0.77):
          * p6 (p = 0.77)
          * p7 (p = 0.77)
          * p8 (p = 0.77)
          * p9 (p = 0.77)
      
      --------------------------------------------------------------------------------
      Note: All tests follow the intersection-union principle and the respective overall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      summary(nine_res_v1)
    Output
      -------------------------------------------------------------------------------- 
      
      Summary of functional redundancy test for p1, p2, p3, p4, and p5.
      
      -------------------------------------------------------------------------------- 
      
      Hypotheses:
      
        H0: The species are not functionally redundant.
        HA: The species are functionally redundant.
      
      Test specification:
      
        Alpha: 0.05
        Confidence level: 90%
        Equivalence margin: (-2, 2) [user specified delta]
      
      Overall result:
      
        Functional redundancy not established (p = 1.00).
      
      Individual Criteria:
      
        Monoculture performances (p = 1.00): [1/10] contrast equivalent to 0.
      
          Contrast  Estimate  p.value
          --------  --------  -------
          p1 - p2   1.110     0.12   
          p1 - p3   1.483     0.24   
          p1 - p4   3.385     0.97   
          p1 - p5   -1.153    0.13   
          p2 - p3   0.373     0.015  
          p2 - p4   2.275     0.64   
          p2 - p5   -2.263    0.64   
          p3 - p4   1.901     0.45   
          p3 - p5   -2.636    0.80   
          p4 - p5   -4.537    1.00   
      
        Within-group interactions (p = 0.31): [0/1] contrasts equivalent to 0.
      
          Contrast        Estimate  p.value
          --------------  --------  -------
          equi - mono_av  -1.602    0.31   
      
        Between-group interactions (p = 0.77): [32/40] contrasts equivalent to 0.
          Showing 20 of 40 contrasts. Use `max.print` to show more.
      
          Contrast       Estimate  p.value
          -------------  --------  -------
          p6_p1 - p6_p2  0.555     0.0001 
          p6_p1 - p6_p3  0.742     0.0004 
          p6_p1 - p6_p4  1.692     0.20   
          p6_p1 - p6_p5  -0.576    0.0001 
          p6_p2 - p6_p3  0.187     0.0001 
          p6_p2 - p6_p4  1.137     0.010  
          p6_p2 - p6_p5  -1.131    0.010  
          p6_p3 - p6_p4  0.951     0.0024 
          p6_p3 - p6_p5  -1.318    0.033  
          p6_p4 - p6_p5  -2.269    0.77   
          p7_p1 - p7_p2  0.555     0.0001 
          p7_p1 - p7_p3  0.742     0.0004 
          p7_p1 - p7_p4  1.692     0.20   
          p7_p1 - p7_p5  -0.576    0.0001 
          p7_p2 - p7_p3  0.187     0.0001 
          p7_p2 - p7_p4  1.137     0.010  
          p7_p2 - p7_p5  -1.131    0.010  
          p7_p3 - p7_p4  0.951     0.0024 
          p7_p3 - p7_p5  -1.318    0.033  
          p7_p4 - p7_p5  -2.269    0.77   
          ... 20 additional contrasts not shown.
      
      Note:  All tests follow the intersection-union principle and the respectiveoverall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      print(nine_res_v2)
    Output
      -------------------------------------------------------------------------------- 
      
      Testing functional redundancy for species p1, p2, p3, p4, p5, p6, p7, p8, and p9.
      
      -------------------------------------------------------------------------------- 
      
      Equivalence tests conducted at the alpha = 0.05 level (i.e, using 90% CIs) with a delta margin of (-2, 2) as defined by user.
      
      Null hypothesis (H0): The species are not functionally redundant.
      Alternative hypothesis (H1): The species are functionally redundant.
      
      --------------------------------------------------------------------------------
      
      Results: 
      
      Functional redundancy was not established for p1, p2, p3, p4, p5, p6, p7, p8, and p9 (p = 1.00) due to the following reasons.
      
        - Monoculture performances were not found to be equivalent (p = 1.00).
        - Within-group interactions were not found to be equivalent to 0 (p = 0.46).
      
      --------------------------------------------------------------------------------
      Note: All tests follow the intersection-union principle and the respective overall p-value is the maximum p-value across all lower level constituent tests.

---

    Code
      summary(nine_res_v2)
    Output
      -------------------------------------------------------------------------------- 
      
      Summary of functional redundancy test for p1, p2, p3, p4, p5, p6, p7, p8, and p9.
      
      -------------------------------------------------------------------------------- 
      
      Hypotheses:
      
        H0: The species are not functionally redundant.
        HA: The species are functionally redundant.
      
      Test specification:
      
        Alpha: 0.05
        Confidence level: 90%
        Equivalence margin: (-2, 2) [user specified delta]
      
      Overall result:
      
        Functional redundancy not established (p = 1.00).
      
      Individual Criteria:
      
        Monoculture performances (p = 1.00): [6/36] contrasts equivalent to 0.
          Showing 20 of 36 contrasts. Use `max.print` to show more.
      
          Contrast  Estimate  p.value
          --------  --------  -------
          p1 - p2   1.110     0.12   
          p1 - p3   1.483     0.24   
          p1 - p4   3.385     0.97   
          p1 - p5   -1.153    0.13   
          p1 - p6   4.214     1.00   
          p1 - p7   4.692     1.00   
          p1 - p8   1.528     0.29   
          p1 - p9   0.622     0.050  
          p2 - p3   0.373     0.015  
          p2 - p4   2.275     0.64   
          p2 - p5   -2.263    0.64   
          p2 - p6   3.104     0.91   
          p2 - p7   3.582     0.97   
          p2 - p8   0.418     0.030  
          p2 - p9   -0.488    0.036  
          p3 - p4   1.901     0.45   
          p3 - p5   -2.636    0.80   
          p3 - p6   2.731     0.81   
          p3 - p7   3.209     0.93   
          p3 - p8   0.045     0.0098 
          ... 16 additional contrasts not shown.
      
        Within-group interactions (p = 0.46): [0/1] contrasts equivalent to 0.
      
          Contrast        Estimate  p.value
          --------------  --------  -------
          equi - mono_av  1.955     0.46   
      
      Note:  All tests follow the intersection-union principle and the respectiveoverall p-value is the maximum p-value across all lower level constituent tests.

