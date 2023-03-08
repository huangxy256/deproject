from deproject.ExampleFolder.example_file import example_function


def test_example_function():

    a = 1
    out = example_function(a=a)
    assert out == 1


