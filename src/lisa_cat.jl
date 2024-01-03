using Catlab.CategoricalAlgebra

@present SchGraph(FreeSchema) begin
  V::Ob
  E::Ob
  src::Hom(E,V)
  tgt::Hom(E,V)
end

@acset_type Graph(SchGraph, index=[:src,:tgt])