Introduction
============

`Eawag envipath_reaction`
-----------------------------

envipath_reaction is a side project of enviRule (https://github.com/zhangky12/enviRule) to test automatically extracted rules on targeted reactions. 
"enviRule: An End-to-end System for Automatic Extraction of Reaction Patterns from Environmental Contaminant Biotransformation Pathways" (https://doi.org/10.1093/bioinformatics/btad407) by Kunyang Zhang and Prof. Dr. Kathrin Fenner.


Note: you can find reaction groups and rules to play with under the folder "example"

`Quick start`
-----------------------------

Make sure you have java installed on the computer, for example:

```
java version "16.0.1" 2021-04-20
Java(TM) SE Runtime Environment (build 16.0.1+9-24)
Java HotSpot(TM) 64-Bit Server VM (build 16.0.1+9-24, mixed mode, sharing)
```

Please always replace "<...>" with your local address in the following examples. The envipath_reaction.jar can be found under out/artifacts/envipath_reaction_jar/

```
java -jar envipath_reaction.jar -reactions <path to any clustered group. E.g., 1-2.txt> -rules <path to the corresponding rule file. E.g., rule-1-2.txt>
```

The rule from your selected rule file will be triggered on each substrate and you will see the products, then you can compare them with the true products in reactions. For example:

In 1-2.txt, we have two reactions:

```
CC1=NC(=C(CN)C=N1)N>>CC1=NC(=C(C=N1)CO)N
C(=C\C(=O)[O-])/C(=C\C(=O)C(=O)[O-])/N>>C(=C\C(=O)[O-])/C(=C\C(=O)C(=O)[O-])/O
```

The products after applying rules are:

```
Reactant: 
CC1=NC(=C(CN)C=N1)N
Products:
CC1=NC(=C(C=N1)CO)N

Reactant: 
C(=C\C(=O)[O-])/C(=C\C(=O)C(=O)[O-])/N
Products:
C(=CC(=O)[O-])C(=O)CC(=O)C(=O)[O-]
```

ps: The second product after applying rules is actually equivalent with the product in the second reaction (enol -> keto). Because products are normalized in envipath_reaction to get rid of ambiguity. Normalized rules can be found in src/main/java/reactor.java, including basicRuleSet, enhancedRuleSet, exoticRuleSet, cutCoARuleSet, and enolKetoRuleSet. 
