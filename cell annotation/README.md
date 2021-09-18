#Cell type annotation

随着单细胞技术的发展，获得单细胞数据越来越成为可能，随着单细胞数据的增多，单细胞亚群的鉴定越来越成为一个单细胞
分析的一个关键的，限制性的过程。面对大量的单细胞数据，就犹如面对天书，不知所云，应用就变得更不可能，因此单细胞
鉴定是把这些数据分类，通过分类鉴定成我们熟悉的细胞类型，我们才能了解细胞的功能和特性，才能解读这种单细胞天书

##细胞鉴定的原因
细胞异质性： 不同的细胞具有不同的功能，不同的表达模式或者说不同的基因组表观
细胞的共性： 每个细胞都有独特的表达模式，都有自己特有的基因，同一类型的细胞都有类似的表达模式
单细胞数据库： 基于已有大量的表达数据库，可以根据基因特征来判断细胞类型，知道什么样的细胞对应什么样的基因表达

##单细胞亚群鉴定方法
- 经典：通过查看亚群特异性表达基因，根据marker基因数据库确定细胞亚群
- 相似性（2019年以后）：通过比较单细胞基因表达与纯特定的细胞类型的bulk rna进行相关性分析，相关性越高，表明细胞类型越相似，比如singleR/scMCA/scHCA属于此类。Seurat也可以做类似事情。
- 半监督： 这是通过给定数据集进行训练，然后根据训练的模型预测新的单细胞数据的细胞类型， 比如cellassign/scibet都属于此类
- 差异基因：这里主要是通过差异基因对marker数据进行富集分析，根据富集程度判断细胞类型，比如clusterprofiler和scsa属于此类
通过使用不同的方法或者集中方法鉴定细胞亚群，确认细胞类型，从而确定细胞的功能特性。

##常见单细胞亚群自动化鉴定的方法

###scHCA/scMCA
通过人的单细胞图谱进行细胞类型自动鉴定，目前有单细胞图谱的物种主要是人和小鼠两个物种，图谱数据既可以当成marker基因数据库来用，也可以直接进行细胞相似性亚群鉴定分析。速度较快。
###singleR v1.0.6
通过比较不同数据集或者bulk RNA的相似性，进行单细胞亚群鉴定。速度快，自动化。
###cellassign通过给定细胞类型的marker基因，通过TensorFlow机器学习的方法，给每种细胞类型进行打分，从而进行细胞类型鉴定。
###clusterprofiler通过使用亚群特异性的高表达的基因进行基因集富集，通过富集程度进行细胞亚群鉴定

##单细胞marker基因数据库
###传统蛋白marker基因
传统蛋白的marker手册，比如BD公司的Human and Mouse CD Marker Handbook和Invitrogen公司的Immune cell guide都是非常不错的手册，免疫细胞几乎都可以使用以上marker基因亚群鉴定
###Cellmarker数据库 
13605 cell markers of 467 cell types in 158 human tissue/sub-tissues and 9148 cell markers of 389 cell types in 81 mouse tissues/subtissue
###panglaodb
8286 association （178 cell types， 4679 gene symbols， 29 tissues）
###单细胞marker基因数据库
- 单细胞图谱（人和小鼠）
https://db.cngb.org/HCL/#  http://bis.zju.edu.cn/MCA/blast.html 小鼠50多种组织的四十万个细胞的单细胞转录组数据，鉴定除了小鼠的98个主要细胞亚群，和800多细胞亚型。
2020年，702968个单细胞的转录组数据，研究团队系统建立了跨越人类胚胎和成年两个时期的单细胞图谱，鉴定了102种细胞簇和843种细胞亚类。
单细胞marker基因数据库-EMBL-EBI https://www.ebi.ac.uk/gxa/sc/home 
单细胞marker基因数据库-Human Cell Atlas https://data.humancellatlas.org/

##singleR
singleR V1.0.0之前的版本，只需要选择human或者mouse两个物种进行分析即可，无需设置参考数据集，参考数据集都是内置在包中，但是运行极慢，及其占用资源；
V1.0.0版本参考数据集可以自己选择，参考数据集更友好，运行速度极快。
SingleR的参考数据集： HumanPrimaryCellAtlasData, BlueprintEncodeData, DatabaseimmuneCellExpressionData(), NovershtemHematopoleticData, MonocoImmuneData, ImmGeneData, MouseRNAseqData
上述7个数据集，人有5个，小鼠有两个，这些数据集有不同的特征，如果是做免疫相关的研究，则选择免疫相关的数据集，细胞亚群鉴定的结果可能会更细致，比如可以鉴定到CD4 naive 或者effector细胞亚群。

##clusterprofiler
clusterprofiler是Y叔开发的目前功能富集最为广泛的工具，由于有marker基因数据库，因此可以用marker基因数据库对单细胞转录组数据的marker gene list进行细胞类型富集分析，根据富集结果进行细胞亚群的判断
#这是做细胞亚群分析最快的方法，但是此方法有个缺陷，那就是如果经过sort后的细胞或者其他纯化的细胞，关键基因不在marker列表中，其结果将极其的不准确。


##单细胞亚群鉴定工具的选择
人和小鼠，如果是人和小鼠的话，建议可以优先考虑使用scMCA和scHCL两种工具进行自动化鉴定，毕竟细胞类型，组织类型较全。
有参考数据集的，可以考虑使用singleR进行鉴定，singleR对于参考数据集的处理方式较为简单，而且运行速度也较快。
有marker gene， 如果有marker基因list，那么建议可以使用cellassign工具进行单细胞类型自动化鉴定。如果有marker基因，而且还有正负表达或者表达阈值的话，可以选择Garnett
其他物种，如果是非人和小鼠的物种，还是需要手动通过经典的鉴定方式进行细胞亚群鉴定。

##单细胞亚群鉴定工具的发展
工具自动化，目前已经有了30多种自动化单细胞自动鉴定工具，因此就工具数目来说，目前已经足够多，但是有的工具不够自动化，今后自动化将会越来越高。
准确性，目前单细胞自动化鉴定工具还存在准确性有待提高的问题，不同的工具鉴定结果差别较大，因此今后的单细胞自动化鉴定工具的准确性需要较大的提高。
效率，目前单细胞自动化鉴定工具准确性较高的工具，其性能有待提高，比如cellassign通过TensorFlow机器学习需要消耗大量的资源和时间。
适应性，目前工具大多数都是使用的人和小鼠两个物种进行的单细胞自动化鉴定，将来学要对其他物种进行覆盖
界面化，界面化的工具进行单细胞鉴定，降低单细胞分析门槛，虽然目前也有界面化的工具，但是大多较简单，数据量过大或者复杂就不好操作。