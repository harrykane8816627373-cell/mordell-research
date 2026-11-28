import Mordell.InnerCoverData

/-! Generated independent kernel decisions for the signed cover. -/

namespace Mordell.InnerCoverData

set_option maxRecDepth 20000
set_option maxHeartbeats 0

private theorem checked38 : tree38.accept 24 degree (1 / 3) true box38 grid38 = true := by
  decide +kernel

private theorem checked37 : tree37.accept 24 degree (1 / 3) true box37 grid37 = true := by
  decide +kernel

private theorem checked36 : tree36.accept 24 degree (1 / 3) true box36 grid36 = true := by
  decide +kernel

private theorem checked35 : tree35.accept 24 degree (1 / 3) true box35 grid35 = true := by
  have hl : gridSplitOK degree 1 false grid35 grid36 := by decide +kernel
  have hr : gridSplitOK degree 1 true grid35 grid37 := by decide +kernel
  exact GridTree.accept_split hl checked36 hr checked37

private theorem checked34 : tree34.accept 24 degree (1 / 3) true box34 grid34 = true := by
  decide +kernel

private theorem checked33 : tree33.accept 24 degree (1 / 3) true box33 grid33 = true := by
  have hl : gridSplitOK degree 0 false grid33 grid34 := by decide +kernel
  have hr : gridSplitOK degree 0 true grid33 grid35 := by decide +kernel
  exact GridTree.accept_split hl checked34 hr checked35

private theorem checked32 : tree32.accept 24 degree (1 / 3) true box32 grid32 = true := by
  decide +kernel

private theorem checked31 : tree31.accept 24 degree (1 / 3) true box31 grid31 = true := by
  decide +kernel

private theorem checked30 : tree30.accept 24 degree (1 / 3) true box30 grid30 = true := by
  decide +kernel

private theorem checked29 : tree29.accept 24 degree (1 / 3) true box29 grid29 = true := by
  have hl : gridSplitOK degree 2 false grid29 grid30 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid29 grid31 := by decide +kernel
  exact GridTree.accept_split hl checked30 hr checked31

private theorem checked28 : tree28.accept 24 degree (1 / 3) true box28 grid28 = true := by
  decide +kernel

private theorem checked27 : tree27.accept 24 degree (1 / 3) true box27 grid27 = true := by
  have hl : gridSplitOK degree 0 false grid27 grid28 := by decide +kernel
  have hr : gridSplitOK degree 0 true grid27 grid29 := by decide +kernel
  exact GridTree.accept_split hl checked28 hr checked29

private theorem checked26 : tree26.accept 24 degree (1 / 3) true box26 grid26 = true := by
  decide +kernel

private theorem checked25 : tree25.accept 24 degree (1 / 3) true box25 grid25 = true := by
  have hl : gridSplitOK degree 1 false grid25 grid26 := by decide +kernel
  have hr : gridSplitOK degree 1 true grid25 grid27 := by decide +kernel
  exact GridTree.accept_split hl checked26 hr checked27

private theorem checked24 : tree24.accept 24 degree (1 / 3) true box24 grid24 = true := by
  decide +kernel

private theorem checked23 : tree23.accept 24 degree (1 / 3) true box23 grid23 = true := by
  decide +kernel

private theorem checked22 : tree22.accept 24 degree (1 / 3) true box22 grid22 = true := by
  decide +kernel

private theorem checked21 : tree21.accept 24 degree (1 / 3) true box21 grid21 = true := by
  have hl : gridSplitOK degree 1 false grid21 grid22 := by decide +kernel
  have hr : gridSplitOK degree 1 true grid21 grid23 := by decide +kernel
  exact GridTree.accept_split hl checked22 hr checked23

private theorem checked20 : tree20.accept 24 degree (1 / 3) true box20 grid20 = true := by
  have hl : gridSplitOK degree 2 false grid20 grid21 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid20 grid24 := by decide +kernel
  exact GridTree.accept_split hl checked21 hr checked24

private theorem checked19 : tree19.accept 24 degree (1 / 3) true box19 grid19 = true := by
  decide +kernel

private theorem checked18 : tree18.accept 24 degree (1 / 3) true box18 grid18 = true := by
  decide +kernel

private theorem checked17 : tree17.accept 24 degree (1 / 3) true box17 grid17 = true := by
  have hl : gridSplitOK degree 2 false grid17 grid18 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid17 grid19 := by decide +kernel
  exact GridTree.accept_split hl checked18 hr checked19

private theorem checked16 : tree16.accept 24 degree (1 / 3) true box16 grid16 = true := by
  have hl : gridSplitOK degree 0 false grid16 grid17 := by decide +kernel
  have hr : gridSplitOK degree 0 true grid16 grid20 := by decide +kernel
  exact GridTree.accept_split hl checked17 hr checked20

private theorem checked15 : tree15.accept 24 degree (1 / 3) true box15 grid15 = true := by
  decide +kernel

private theorem checked14 : tree14.accept 24 degree (1 / 3) true box14 grid14 = true := by
  decide +kernel

private theorem checked13 : tree13.accept 24 degree (1 / 3) true box13 grid13 = true := by
  have hl : gridSplitOK degree 2 false grid13 grid14 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid13 grid15 := by decide +kernel
  exact GridTree.accept_split hl checked14 hr checked15

private theorem checked12 : tree12.accept 24 degree (1 / 3) true box12 grid12 = true := by
  decide +kernel

private theorem checked11 : tree11.accept 24 degree (1 / 3) true box11 grid11 = true := by
  have hl : gridSplitOK degree 0 false grid11 grid12 := by decide +kernel
  have hr : gridSplitOK degree 0 true grid11 grid13 := by decide +kernel
  exact GridTree.accept_split hl checked12 hr checked13

private theorem checked10 : tree10.accept 24 degree (1 / 3) true box10 grid10 = true := by
  have hl : gridSplitOK degree 1 false grid10 grid11 := by decide +kernel
  have hr : gridSplitOK degree 1 true grid10 grid16 := by decide +kernel
  exact GridTree.accept_split hl checked11 hr checked16

private theorem checked9 : tree9.accept 24 degree (1 / 3) true box9 grid9 = true := by
  have hl : gridSplitOK degree 2 false grid9 grid10 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid9 grid25 := by decide +kernel
  exact GridTree.accept_split hl checked10 hr checked25

private theorem checked8 : tree8.accept 24 degree (1 / 3) true box8 grid8 = true := by
  have hl : gridSplitOK degree 2 false grid8 grid9 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid8 grid32 := by decide +kernel
  exact GridTree.accept_split hl checked9 hr checked32

private theorem checked7 : tree7.accept 24 degree (1 / 3) true box7 grid7 = true := by
  decide +kernel

private theorem checked6 : tree6.accept 24 degree (1 / 3) true box6 grid6 = true := by
  have hl : gridSplitOK degree 1 false grid6 grid7 := by decide +kernel
  have hr : gridSplitOK degree 1 true grid6 grid8 := by decide +kernel
  exact GridTree.accept_split hl checked7 hr checked8

private theorem checked5 : tree5.accept 24 degree (1 / 3) true box5 grid5 = true := by
  decide +kernel

private theorem checked4 : tree4.accept 24 degree (1 / 3) true box4 grid4 = true := by
  decide +kernel

private theorem checked3 : tree3.accept 24 degree (1 / 3) true box3 grid3 = true := by
  have hl : gridSplitOK degree 1 false grid3 grid4 := by decide +kernel
  have hr : gridSplitOK degree 1 true grid3 grid5 := by decide +kernel
  exact GridTree.accept_split hl checked4 hr checked5

private theorem checked2 : tree2.accept 24 degree (1 / 3) true box2 grid2 = true := by
  have hl : gridSplitOK degree 0 false grid2 grid3 := by decide +kernel
  have hr : gridSplitOK degree 0 true grid2 grid6 := by decide +kernel
  exact GridTree.accept_split hl checked3 hr checked6

private theorem checked1 : tree1.accept 24 degree (1 / 3) true box1 grid1 = true := by
  have hl : gridSplitOK degree 2 false grid1 grid2 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid1 grid33 := by decide +kernel
  exact GridTree.accept_split hl checked2 hr checked33

private theorem checked0 : tree0.accept 24 degree (1 / 3) true box0 grid0 = true := by
  have hl : gridSplitOK degree 2 false grid0 grid1 := by decide +kernel
  have hr : gridSplitOK degree 2 true grid0 grid38 := by decide +kernel
  exact GridTree.accept_split hl checked1 hr checked38

theorem archived_tree_accepted : tree.accept 24 degree (1 / 3) true root grid0 = true := checked0

end Mordell.InnerCoverData
