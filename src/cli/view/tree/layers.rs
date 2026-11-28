//! Layers drawn over the tree: values as labels, and one value as branch
//! colors. New layers only need a `TreeLayer` from `segul-core` and a slot
//! here; drawing does not change.

use segul_core::core::tree::view::{
    AnnotationLayer, BranchLengthLayer, LayerValue, SupportLayer, TreeLayer, annotation_keys,
};
use segul_core::parser::tree::{NodeId, Tree};

/// A layer and whether its labels are shown.
pub(super) struct LayerSlot {
    pub(super) layer: Box<dyn TreeLayer>,
    pub(super) visible: bool,
}

/// Layers in drawing order, plus the one coloring the branches.
pub(super) struct Layers {
    pub(super) slots: Vec<LayerSlot>,
    /// Slot coloring the branches, with its smallest and largest value.
    color: Option<(usize, f64, f64)>,
}

impl Layers {
    /// Support and branch lengths, then every annotation key in `tree`.
    /// Visibility and branch colors carry over from `previous` by name.
    pub(super) fn for_tree(tree: &Tree, previous: Option<&Layers>) -> Self {
        let mut slots = vec![
            LayerSlot {
                layer: Box::new(SupportLayer),
                visible: true,
            },
            LayerSlot {
                layer: Box::new(BranchLengthLayer),
                visible: false,
            },
        ];
        for key in annotation_keys(tree) {
            slots.push(LayerSlot {
                layer: Box::new(AnnotationLayer::new(&key)),
                visible: false,
            });
        }
        let mut layers = Self { slots, color: None };
        if let Some(previous) = previous {
            for slot in &mut layers.slots {
                if let Some(old) = previous.find(slot.layer.name()) {
                    slot.visible = previous.slots[old].visible;
                }
            }
            if let Some(name) = previous.color_name() {
                let index = layers.find(name);
                layers.set_color(tree, index);
            }
        }
        layers
    }

    fn find(&self, name: &str) -> Option<usize> {
        self.slots.iter().position(|s| s.layer.name() == name)
    }

    /// Show or hide slot `index`.
    pub(super) fn toggle(&mut self, index: usize) -> bool {
        match self.slots.get_mut(index) {
            Some(slot) => {
                slot.visible = !slot.visible;
                true
            }
            None => false,
        }
    }

    /// Name of the layer coloring the branches.
    pub(super) fn color_name(&self) -> Option<&str> {
        self.color.map(|(i, _, _)| self.slots[i].layer.name())
    }

    /// Color branches by the next layer with numbers, then by none.
    pub(super) fn cycle_color(&mut self, tree: &Tree) {
        let start = self.color.map_or(0, |(i, _, _)| i + 1);
        for index in start..self.slots.len() {
            if self.set_color(tree, Some(index)) {
                return;
            }
        }
        self.color = None;
    }

    // Colors by slot `index` when any node has a number for it. Scans the
    // tree once for the value range.
    fn set_color(&mut self, tree: &Tree, index: Option<usize>) -> bool {
        self.color = None;
        let Some(index) = index else {
            return false;
        };
        let layer = &self.slots[index].layer;
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for id in tree.preorder() {
            if let Some(LayerValue::Number(v)) = layer.value(tree, id) {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        if lo > hi {
            return false;
        }
        self.color = Some((index, lo, hi));
        true
    }

    /// Visible values for `id`, joined by slashes.
    pub(super) fn labels(&self, tree: &Tree, id: NodeId) -> Option<String> {
        let values: Vec<String> = self
            .slots
            .iter()
            .filter(|s| s.visible)
            .filter_map(|s| s.layer.value(tree, id))
            .map(|v| v.to_string())
            .collect();
        (!values.is_empty()).then(|| values.join("/"))
    }

    /// Every layer value for `id`, as `name value` pairs.
    pub(super) fn describe(&self, tree: &Tree, id: NodeId) -> Vec<String> {
        self.slots
            .iter()
            .filter_map(|s| {
                s.layer
                    .value(tree, id)
                    .map(|v| format!("{} {}", s.layer.name(), v))
            })
            .collect()
    }

    /// Branch color position of `id`, from 0 to 1.
    pub(super) fn color_at(&self, tree: &Tree, id: NodeId) -> Option<f64> {
        let (index, lo, hi) = self.color?;
        match self.slots[index].layer.value(tree, id)? {
            LayerValue::Number(v) if hi > lo => Some((v - lo) / (hi - lo)),
            LayerValue::Number(_) => Some(1.0),
            LayerValue::Text(_) => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use segul_core::parser::tree::parse_tree_str;

    #[test]
    fn test_layers() {
        let tree = parse_tree_str("((A[&rate=1]:1,B[&rate=3]:2)90:1,C:3);").unwrap();
        let mut layers = Layers::for_tree(&tree, None);
        let names: Vec<&str> = layers.slots.iter().map(|s| s.layer.name()).collect();
        assert_eq!(names, ["Support", "Lengths", "rate"]);

        let a = tree.leaves().next().unwrap();
        let clade = tree.parent(a).unwrap();
        assert_eq!(layers.labels(&tree, clade).as_deref(), Some("90"));
        assert_eq!(layers.labels(&tree, a), None);
        layers.toggle(1);
        assert_eq!(layers.labels(&tree, clade).as_deref(), Some("90/1"));
        assert!(!layers.toggle(9));

        layers.cycle_color(&tree);
        assert_eq!(layers.color_name(), Some("Support"));
        assert_eq!(layers.color_at(&tree, clade), Some(1.0));
        layers.cycle_color(&tree);
        layers.cycle_color(&tree);
        assert_eq!(layers.color_name(), Some("rate"));
        assert_eq!(layers.color_at(&tree, a), Some(0.0));
        layers.cycle_color(&tree);
        assert_eq!(layers.color_name(), None);

        layers.cycle_color(&tree);
        let next = Layers::for_tree(&tree, Some(&layers));
        assert!(next.slots[1].visible);
        assert_eq!(next.color_name(), Some("Support"));
    }
}
