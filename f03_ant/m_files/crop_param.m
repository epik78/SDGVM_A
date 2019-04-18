%Gets the crop parameters to convert production to NPP accrding to crop id
function dis=crop_param(a)

%abaca 
crop_params(1,:)=[.28 .8 .5];
%agave
crop_params(2,:)=[.28 .5 .5];
%alfalfa
crop_params(3,:)=[1 .2 .53];
%almond
crop_params(4,:)=[.28 .9 .75];
%aniseetc
crop_params(5,:)=[.28 .8 .5];
%apple
crop_params(6,:)=[.3 .16 .75];
%apricot
crop_params(7,:)=[.3 .14 .75];
%areca
crop_params(8,:)=[.28 .8 .5];
%artichoke
crop_params(9,:)=[.45 .13 .85];
%asparagus
crop_params(10,:)=[.45 .08 .85];
%avocado
crop_params(11,:)=[.3 .26 .75];
%bambara
crop_params(12,:)=[.49 .9 .85];
%banana
crop_params(13,:)=[.3 .2 .75];
%barley
crop_params(14,:)=[.49 .89 .5];
%bean
crop_params(15,:)=[.45 .1 .85];
%beetfor
crop_params(16,:)=[1 .13 .85];
%berrynes
crop_params(17,:)=[.3 .19 .75];
%blueberry
crop_params(18,:)=[.3 .15 .75];
%brazil
crop_params(19,:)=[.28 .8 .75];
%broadbean
crop_params(20,:)=[.45 .13 .85];
%buckwheat
crop_params(21,:)=[.4 .88 .8];
%cabbage
crop_params(22,:)=[.45 .08 .85];
%cabbagefor
crop_params(23,:)=[1 .08 .85];
%canaryseed
crop_params(24,:)=[.4 .88 .8];
%carob
crop_params(25,:)=[.3 .19 .75];
%carrot
crop_params(26,:)=[.45 .12 .85];
%carrotfor
crop_params(27,:)=[1 .12 .85];
%cashew
crop_params(28,:)=[.28 .8 .75];
%cashewapple
crop_params(29,:)=[.3 .19 .75];
%cassava
crop_params(30,:)=[.48 .32 .85];
%castor
crop_params(31,:)=[.52 .73 .8];
%cauliflower
crop_params(32,:)=[.45 .08 .85];
%cerealnes
crop_params(33,:)=[.4 .88 .8];
%cherry
crop_params(34,:)=[.3 .14 .75];
%chestnut
crop_params(35,:)=[.28 .8 .75];
%chickpea
crop_params(36,:)=[.44 .9 .85];
%chicory
crop_params(37,:)=[.28 .8 .8];
%chilleetc
crop_params(38,:)=[.45 .08 .85];
%cinnamon
crop_params(39,:)=[.28 .8 .5];
%citrusnes
crop_params(40,:)=[.3 .13 .75];
%clove
crop_params(41,:)=[.28 .8 .5];
%clover
crop_params(42,:)=[1 .2 .5];
%cocoa
crop_params(43,:)=[.28 .8 .5];
%coconut
crop_params(44,:)=[.28 .8 .5];
%coffee
crop_params(45,:)=[.28 .8 .5];
%cotton
crop_params(46,:)=[.55 .92 .86];
%cowpea
crop_params(47,:)=[.55 .9 .85];
%cranberry
crop_params(48,:)=[.3 .19 .75];
%cucumberetc
crop_params(49,:)=[.45 .04 .85];
%currant
crop_params(50,:)=[.3 .19 .75];
%date
crop_params(51,:)=[.3 .19 .75];
%eggplant
crop_params(52,:)=[.45 .08 .85];
%fibrenes
crop_params(53,:)=[.28 .8 .8];
%fig
crop_params(54,:)=[.3 .21 .75];
%flax
crop_params(55,:)=[.28 .8 .8];
%fonio
crop_params(56,:)=[.4 .88 .8];
%fornes
crop_params(57,:)=[1 .2 .65];
%fruitnes
crop_params(58,:)=[.3 .19 .75];
%garlic
crop_params(59,:)=[.45 .13 .85];
%ginger
crop_params(60,:)=[.28 .8 .5];
%gooseberry
crop_params(61,:)=[.3 .19 .75];
%grape
crop_params(62,:)=[.3 .19 .75];
%grapefruitetc
crop_params(63,:)=[.3 .09 .75];
%grassnes
crop_params(64,:)=[1 .2 .65];
%greenbean
crop_params(65,:)=[.45 .1 .85];
%greenbroadbean
crop_params(66,:)=[.45 .13 .85];
%greencorn
crop_params(67,:)=[.45 .13 .85];
%greenonion
crop_params(68,:)=[.45 .09 .85];
%greenpea
crop_params(69,:)=[.45 .13 .85];
%groundnut
crop_params(70,:)=[.4 .92 .8];
%hazelnut
crop_params(71,:)=[.28 .8 .75];
%hemp
crop_params(72,:)=[.28 .8 .8];
%hempseed
crop_params(73,:)=[.52 .73 .8];
%hop
crop_params(74,:)=[.28 .8 .5];
%jute
crop_params(75,:)=[.28 .8 .8];
%jutelikefiber
crop_params(76,:)=[.28 .8 .8];
%kapokfiber
crop_params(77,:)=[.28 .8 .5];
%kapokseed
crop_params(78,:)=[.28 .8 .5];
%karite
crop_params(79,:)=[.28 .8 .5];
%kiwi
crop_params(80,:)=[.3 .13 .75];
%kolanut
crop_params(81,:)=[.28 .8 .5];
%legumenes
crop_params(82,:)=[1 .2 .65];
%lemonlime
crop_params(83,:)=[.3 .13 .75];
%lentil
crop_params(84,:)=[.46 .89 .85];
%lettuce
crop_params(85,:)=[.45 .05 .85];
%linseed
crop_params(86,:)=[.52 .73 .8];
%lupin
crop_params(87,:)=[.41 .89 .85];
%maize
crop_params(88,:)=[.45 .89 .85];
%maizefor
crop_params(89,:)=[1 .35 .85];
%mango
crop_params(90,:)=[.3 .19 .75];
%mate
crop_params(91,:)=[.28 .8 .5];
%melonetc
crop_params(92,:)=[.45 .1 .85];
%melonseed
crop_params(93,:)=[.52 .73 .8];
%millet
crop_params(94,:)=[.4 .9 .88];
%mixedgrain
crop_params(95,:)=[.4 .88 .8];
%mixedgrass
crop_params(96,:)=[1 .2 .65];
%mushroom
crop_params(97,:)=[.45 .13 .85];
%mustard
crop_params(98,:)=[.52 .73 .8];
%nutmeg
crop_params(99,:)=[.28 .8 .5];
%nutnes
crop_params(100,:)=[.28 .8 .5];
%oats
crop_params(101,:)=[.4 .89 .71];
%oilpalm
crop_params(102,:)=[.28 .8 .5];
%oilseedfor
crop_params(103,:)=[1 .35 .8];
%oilseednes
crop_params(104,:)=[.52 .73 .8];
%okra
crop_params(105,:)=[.45 .1 .85];
%olive
crop_params(106,:)=[.28 .8 .5];
%onion
crop_params(107,:)=[.45 .13 .85];
%orange
crop_params(108,:)=[.3 .13 .5];
%papaya
crop_params(109,:)=[.3 .11 .75];
%pea
crop_params(110,:)=[.45 .89 .85];
%peachetc
crop_params(111,:)=[.3 .14 .75];
%pear
crop_params(112,:)=[.3 .16 .75];
%pepper
crop_params(113,:)=[.28 .8 .5];
%peppermint
crop_params(114,:)=[.28 .8 .5];
%persimmon
crop_params(115,:)=[.3 .36 .75];
%pigeonpea
crop_params(116,:)=[.23 .9 .85];
%pimento
crop_params(117,:)=[.28 .8 .8];
%pineapple
crop_params(118,:)=[.3 .14 .75];
%pistachio
crop_params(119,:)=[.28 .8 .75];
%plantain
crop_params(120,:)=[.3 .2 .75];
%plum
crop_params(121,:)=[.3 .15 .75];
%poppy
crop_params(122,:)=[.52 .73 .8];
%potato
crop_params(123,:)=[.5 .28 .8];
%pulsenes
crop_params(124,:)=[.49 .9 .85];
%pumpkinetc
crop_params(125,:)=[.45 .2 .85];
%pyrethrum
crop_params(126,:)=[.28 .8 .5];
%quince
crop_params(127,:)=[.3 .16 .75];
%quinoa
crop_params(128,:)=[.4 .88 .8];
%ramie
crop_params(129,:)=[.28 .8 .5];
%rapeseed
crop_params(130,:)=[.3 .73 .8];
%rasberry
crop_params(131,:)=[.3 .13 .75];
%rice
crop_params(132,:)=[.4 .89 .8];
%rootnes
crop_params(133,:)=[.4 .2 .8];
%rubber
crop_params(134,:)=[.28 .8 .5];
%rye
crop_params(135,:)=[.35 .88 .76];
%ryefor
crop_params(136,:)=[1 .2 .65];
%safflower
crop_params(137,:)=[.52 .91 .8];
%sesame
crop_params(138,:)=[.52 .92 .8];
%sisal
crop_params(139,:)=[.28 .8 .5];
%sorghum
crop_params(140,:)=[.4 .89 .8];
%sorghumfor
crop_params(141,:)=[1 .35 .85];
%sourcherry
crop_params(142,:)=[.3 .14 .75];
%soybean
crop_params(143,:)=[.42 .91 .85];
%spicenes
crop_params(144,:)=[.28 .8 .5];
%spinach
crop_params(145,:)=[.45 .08 .85];
%stonefruitnes
crop_params(146,:)=[.3 .19 .75];
%strawberry
crop_params(147,:)=[.3 .08 .75];
%stringbean
crop_params(148,:)=[.45 .13 .85];
%sugarbeet
crop_params(149,:)=[.4 .12 .8];
%sugarcane
crop_params(150,:)=[.85 .15 .85];
%sugarnes
crop_params(151,:)=[.28 .56 .85];
%sunflower
crop_params(152,:)=[.39 .94 .94];
%swedefor
crop_params(153,:)=[1 .13 .85];
%sweetpotato
crop_params(154,:)=[.5 .25 .8];
%tangetc
crop_params(155,:)=[.3 .19 .75];
%taro
crop_params(156,:)=[.4 .2 .8];
%tea
crop_params(157,:)=[.28 .8 .5];
%tobacco
crop_params(158,:)=[.28 .8 .8];
%tomato
crop_params(159,:)=[.45 .06 .85];
%triticale
crop_params(160,:)=[.46 .9 .8];
%tropicalnes
crop_params(161,:)=[.3 .19 .75];
%tung
crop_params(162,:)=[.28 .8 .5];
%turnipfor
crop_params(163,:)=[1 .13 .85];
%vanilla
crop_params(164,:)=[.28 .8 .5];
%vegetablenes
crop_params(165,:)=[.45 .13 .85];
%vegfor
crop_params(166,:)=[1 .13 .85];
%vetch
crop_params(167,:)=[.49 .9 .85];
%walnut
crop_params(168,:)=[.28 .91 .75];
%watermelon
crop_params(169,:)=[.45 .09 .85];
%wheat
crop_params(170,:)=[.39 .89 .81];
%yam
crop_params(171,:)=[.4 .3 .8];
%yautia
crop_params(172,:)=[.4 .2 .8];

dis=crop_params(a,:);    


end
