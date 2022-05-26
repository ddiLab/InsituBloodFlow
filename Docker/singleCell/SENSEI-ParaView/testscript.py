import paraview.simple


def main():
    render_view = paraview.simple.GetActiveViewOrCreate('RenderView')
    render_view.ViewSize = [1920, 1080]

    sphere = paraview.simple.Sphere()
    sphere_display = paraview.simple.Show(sphere, render_view)
    sphere_display.DiffuseColor = [1.0, 0.0, 0.0]

    render_view.CameraPosition = [-5.0, 0.0, 0.0]
    render_view.CameraFocalPoint = [0.0, 0.0, 0.0]
    render_view.CameraViewUp = [0.0, 0.0, 1.0]
    render_view.CameraParallelScale = 1.0

    paraview.simple.SaveScreenshot('./test.png', render_view)


if __name__ == '__main__':
    main()
