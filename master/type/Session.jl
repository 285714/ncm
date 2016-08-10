"""
Comprises everything needed for path following: `Project`, `SystemCore`,
`ContinuationMethod`, and `Visualization`. Sessions can be managed via `create`, `save`,
and `load`.
"""
type Session
	P::Project
	core::SystemCore
	cont::ContinuationMethod
	viz::Visualization
	Session() = new()
end
